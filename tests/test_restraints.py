# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import pytest
import torch

from chai_lab import chai1
from chai_lab.data.collate.collate import Collate
from chai_lab.data.dataset.all_atom_feature_context import AllAtomFeatureContext
from chai_lab.data.dataset.constraints.restraint_context import (
    load_manual_restraints_for_chai1,
)
from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.inference_dataset import Input, load_chains_from_raw
from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.templates.context import TemplateContext
from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.data.parsing.restraints import (
    PairwiseInteraction,
    PairwiseInteractionType,
    parse_pairwise_table,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.paths import repo_root


def test_loading_restraints():
    """Small test to ensure that restraints can be loaded."""
    contact_path = repo_root / "examples" / "restraints" / "contact.restraints"
    pocket_path = repo_root / "examples" / "restraints" / "pocket.restraints"

    assert len(parse_pairwise_table(contact_path)) > 0
    assert len(parse_pairwise_table(pocket_path)) > 0


@pytest.mark.parametrize(
    "entity_name_as_subchain,restraints_wrt_entity_names",
    [
        (True, True),  # Both given w.r.t. automatic names; should load correctly
        (True, False),  # Mismatched; should not load
        (False, True),  # Mismatched; should not load
        (False, False),  # Both w.r.t. fasta-derived entity names; should load correctly
    ],
)
def test_restraints_with_manual_chain_names(
    entity_name_as_subchain: bool, restraints_wrt_entity_names: bool
):
    """subchain ID scheme and restraint scheme compatibility"""
    inputs = [
        Input("GGGGGG", entity_type=EntityType.PROTEIN.value, entity_name="G"),
        Input("HHHHHH", entity_type=EntityType.PROTEIN.value, entity_name="H"),
    ]

    restraints = [
        PairwiseInteraction(
            chainA="G" if restraints_wrt_entity_names else "A",
            res_idxA="G1",
            atom_nameA="",
            chainB="H" if restraints_wrt_entity_names else "B",
            res_idxB="H1",
            atom_nameB="",
            connection_type=PairwiseInteractionType.CONTACT,
        ),
        PairwiseInteraction(
            chainA="G" if restraints_wrt_entity_names else "A",
            res_idxA="",
            atom_nameA="",
            chainB="H" if restraints_wrt_entity_names else "B",
            res_idxB="H1",
            atom_nameB="",
            connection_type=PairwiseInteractionType.POCKET,
        ),
    ]

    chains = load_chains_from_raw(
        inputs=inputs, entity_name_as_subchain=entity_name_as_subchain
    )
    assert len(chains) == 2

    structure_context = AllAtomStructureContext.merge(
        [c.structure_context for c in chains]
    )
    ft_ctx = AllAtomFeatureContext(
        chains=chains,
        structure_context=structure_context,
        msa_context=MSAContext.create_single_seq(
            dataset_source=MSADataSource.QUERY,
            tokens=structure_context.token_residue_type.to(dtype=torch.uint8),
        ),
        profile_msa_context=MSAContext.create_single_seq(
            dataset_source=MSADataSource.QUERY,
            tokens=structure_context.token_residue_type.to(dtype=torch.uint8),
        ),
        template_context=TemplateContext.empty(
            n_templates=1, n_tokens=structure_context.num_tokens
        ),
        embedding_context=EmbeddingContext.empty(n_tokens=structure_context.num_tokens),
        restraint_context=load_manual_restraints_for_chai1(chains, None, restraints),
    )

    collator = Collate(
        feature_factory=chai1.feature_factory, num_key_atoms=128, num_query_atoms=32
    )

    batch = collator([ft_ctx])

    assert batch
    ft = batch["features"]
    contact_ft = ft["TokenDistanceRestraint"]
    contact_ft_all_null = torch.allclose(contact_ft, torch.tensor(-1).float())
    pocket_ft = ft["TokenPairPocketRestraint"]
    pocket_ft_all_null = torch.allclose(pocket_ft, torch.tensor(-1).float())

    if entity_name_as_subchain == restraints_wrt_entity_names:
        # Loaded correctly, so some should not be null
        assert not contact_ft_all_null
        assert not pocket_ft_all_null
    else:
        # Did not load; all null
        assert contact_ft_all_null
        assert pocket_ft_all_null

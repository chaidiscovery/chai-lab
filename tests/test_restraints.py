# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging

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


def _restraint_features(
    entity_name_as_subchain: bool, restraints_wrt_entity_names: bool
):
    """Build the collated features for a two-chain input with manual restraints.

    When the two flags disagree the restraints refer to chain names that do not
    exist in the input, so they cannot be placed.
    """
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
    return batch["features"]


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
    ft = _restraint_features(entity_name_as_subchain, restraints_wrt_entity_names)
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


def test_unplaceable_restraints_say_the_run_is_unrestrained(caplog):
    """Restraints that cannot be placed must say so, and say what it costs.

    The features themselves fall back to the null value, which is what
    test_restraints_with_manual_chain_names pins and is the intended behavior.
    The problem this covers is that the null feature is indistinguishable from
    the one produced when no restraints were requested at all, so the log is the
    only signal a caller has that the restraints they asked for were dropped.
    """
    with caplog.at_level(logging.DEBUG):
        ft = _restraint_features(
            entity_name_as_subchain=True, restraints_wrt_entity_names=False
        )

    # Precondition: these restraints really did fail to land.
    assert torch.allclose(ft["TokenPairPocketRestraint"], torch.tensor(-1).float())
    assert torch.allclose(ft["TokenDistanceRestraint"], torch.tensor(-1).float())

    failures = [r for r in caplog.records if r.levelno >= logging.ERROR]
    assert failures, "dropping the requested restraints was not reported at all"
    # The consequence has to be in the message, not left for the reader to infer
    # from a bare exception string.
    assert any("NO pocket restraints" in r.getMessage() for r in failures)
    assert any("NO distance restraints" in r.getMessage() for r in failures)
    # logger.exception keeps the traceback, so the failing assert is identifiable.
    assert any(r.exc_info is not None for r in failures)
    # And nothing may claim success.
    assert not any("Generated pocket feature" in r.getMessage() for r in caplog.records)


def test_placeable_restraints_report_success(caplog):
    """The success line fires only after the restraints are actually applied."""
    with caplog.at_level(logging.DEBUG):
        ft = _restraint_features(
            entity_name_as_subchain=True, restraints_wrt_entity_names=True
        )

    assert not torch.allclose(ft["TokenPairPocketRestraint"], torch.tensor(-1).float())
    assert not torch.allclose(ft["TokenDistanceRestraint"], torch.tensor(-1).float())
    assert any("Generated pocket feature" in r.getMessage() for r in caplog.records)
    assert any("Generated distance feature" in r.getMessage() for r in caplog.records)
    assert not [r for r in caplog.records if r.levelno >= logging.ERROR]

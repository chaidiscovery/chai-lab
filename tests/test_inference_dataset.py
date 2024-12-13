# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Tests for inference dataset.
"""

import pytest
import torch

from chai_lab.data.dataset.inference_dataset import Input, load_chains_from_raw
from chai_lab.data.dataset.structure.all_atom_residue_tokenizer import (
    AllAtomResidueTokenizer,
)
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.sources.rdkit import RefConformerGenerator


@pytest.fixture
def tokenizer() -> AllAtomResidueTokenizer:
    return AllAtomResidueTokenizer(RefConformerGenerator())


def test_malformed_smiles(tokenizer: AllAtomResidueTokenizer):
    """Malformed SMILES should be dropped."""
    # Zn ligand is malformed (should be [Zn+2])
    inputs = [
        Input("RKDESES", entity_type=EntityType.PROTEIN.value, entity_name="foo"),
        Input("Zn", entity_type=EntityType.LIGAND.value, entity_name="bar"),
        Input("RKEEE", entity_type=EntityType.PROTEIN.value, entity_name="baz"),
        Input("EEEEEEEEEEEE", entity_type=EntityType.PROTEIN.value, entity_name="boz"),
    ]
    chains = load_chains_from_raw(
        inputs,
        identifier="test",
        tokenizer=tokenizer,
    )
    assert len(chains) == 3
    for chain in chains:
        # NOTE this check is only valid because there are no residues that are tokenized per-atom
        # Ensures that the entity data and the structure context in each chain are paired correctly
        assert chain.structure_context.num_tokens == len(
            chain.entity_data.full_sequence
        )


def test_ions_parsing(tokenizer: AllAtomResidueTokenizer):
    """Ions as SMILES strings should carry the correct charge."""
    inputs = [Input("[Mg+2]", entity_type=EntityType.LIGAND.value, entity_name="foo")]
    chains = load_chains_from_raw(inputs, identifier="foo", tokenizer=tokenizer)
    assert len(chains) == 1
    chain = chains[0]
    assert chain.structure_context.num_atoms == 1
    assert chain.structure_context.atom_ref_charge == 2
    assert chain.structure_context.atom_ref_element.item() == 12


def test_protein_with_smiles(tokenizer: AllAtomResidueTokenizer):
    """Complex with multiple duplicated protein chains and SMILES ligands."""
    # Based on https://www.rcsb.org/structure/1AFS
    seq = "MDSISLRVALNDGNFIPVLGFGTTVPEKVAKDEVIKATKIAIDNGFRHFDSAYLYEVEEEVGQAIRSKIEDGTVKREDIFYTSKLWSTFHRPELVRTCLEKTLKSTQLDYVDLYIIHFPMALQPGDIFFPRDEHGKLLFETVDICDTWEAMEKCKDAGLAKSIGVSNFNCRQLERILNKPGLKYKPVCNQVECHLYLNQSKMLDYCKSKDIILVSYCTLGSSRDKTWVDQKSPVLLDDPVLCAIAKKYKQTPALVALRYQLQRGVVPLIRSFNAKRIKELTQVFEFQLASEDMKALDGLNRNFRYNNAKYFDDHPNHPFTDEN"
    nap = "NC(=O)c1ccc[n+](c1)[CH]2O[CH](CO[P]([O-])(=O)O[P](O)(=O)OC[CH]3O[CH]([CH](O[P](O)(O)=O)[CH]3O)n4cnc5c(N)ncnc45)[CH](O)[CH]2O"
    tes = "O=C4C=C3C(C2CCC1(C(CCC1O)C2CC3)C)(C)CC4"
    inputs = [
        Input(seq, EntityType.PROTEIN.value, entity_name="A"),
        Input(seq, EntityType.PROTEIN.value, entity_name="B"),
        Input(nap, EntityType.LIGAND.value, entity_name="C"),
        Input(nap, EntityType.LIGAND.value, entity_name="D"),
        Input(tes, EntityType.LIGAND.value, entity_name="E"),
        Input(tes, EntityType.LIGAND.value, entity_name="F"),
    ]
    chains: list[Chain] = load_chains_from_raw(inputs, tokenizer=tokenizer)
    assert len(chains) == len(inputs)

    example = AllAtomStructureContext.merge(
        [chain.structure_context for chain in chains]
    )

    # Should be 1 protein chain, 2 ligand chains
    assert example.token_entity_id.unique().numel() == 3
    assert example.token_asym_id.unique().numel() == 6

    # Check protein chains
    prot_entity_ids = example.token_entity_id[
        example.token_entity_type == EntityType.PROTEIN.value
    ]
    assert torch.unique(prot_entity_ids).numel() == 1
    prot_sym_ids = example.token_sym_id[
        example.token_entity_type == EntityType.PROTEIN.value
    ]
    assert torch.unique(prot_sym_ids).numel() == 2  # Two copies of this chain

    # Check ligand chains
    lig_entity_ids = example.token_entity_id[
        example.token_entity_type == EntityType.LIGAND.value
    ]
    assert torch.unique(lig_entity_ids).numel() == 2
    lig_sym_ids = example.token_sym_id[
        example.token_entity_type == EntityType.LIGAND.value
    ]
    assert torch.unique(lig_sym_ids).numel() == 2  # Two copies of each ligand

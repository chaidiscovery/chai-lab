# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Tests for inference dataset.
"""

from chai_lab.data.dataset.inference_dataset import Input, load_chains_from_raw
from chai_lab.data.dataset.structure.all_atom_residue_tokenizer import (
    AllAtomResidueTokenizer,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.sources.rdkit import RefConformerGenerator


def test_malformed_smiles():
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
        tokenizer=AllAtomResidueTokenizer(RefConformerGenerator()),
    )
    assert len(chains) == 3
    for chain in chains:
        # NOTE this check is only valid because there are no residues that are tokenized per-atom
        # Ensures that the entity data and the structure context in each chain are paired correctly
        assert chain.structure_context.num_tokens == len(
            chain.entity_data.full_sequence
        )

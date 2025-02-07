# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from chai_lab.data.sources.rdkit import RefConformerGenerator


def test_ref_conformer_from_smiles():
    """Test ref conformer generation from SMILES."""
    smiles = "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C"
    rcg = RefConformerGenerator()

    conformer = rcg.generate(smiles)

    assert len(set(conformer.atom_names)) == conformer.num_atoms


def test_ref_conformer_glycan_ccd():
    """Ref conformer from CCD code for a sugar ring."""
    rcg = RefConformerGenerator()
    conformer = rcg.get("MAN")
    assert conformer is not None

    assert len(set(conformer.atom_names)) == conformer.num_atoms

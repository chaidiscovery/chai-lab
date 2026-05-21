# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from chai_lab.chai1 import make_all_atom_feature_context
from chai_lab.data.parsing.glycans import _glycan_string_to_sugars_and_bonds


@pytest.mark.parametrize("ccd_code", ["MAN", "99K", "FUC"])
def test_parsing_ccd_codes(ccd_code: str):
    """Test that various single CCD codes are parsed correctly."""
    res, _ = _glycan_string_to_sugars_and_bonds(ccd_code)
    assert len(res) == 1


def test_complex_parsing():
    glycan = "MAN(6-1 FUC)(4-1 MAN(6-1 MAN(6-1 MAN)))".replace(" ", "")
    sugars, bonds = _glycan_string_to_sugars_and_bonds(glycan)
    assert len(sugars) == 5

    bond1, bond2, bond3, bond4 = bonds

    assert bond1.src_sugar_index == 0
    assert bond1.dst_sugar_index == 1
    assert bond1.src_atom == 6
    assert bond1.dst_atom == 1
    assert bond2.src_sugar_index == 0
    assert bond2.dst_sugar_index == 2
    assert bond2.src_atom == 4
    assert bond2.dst_atom == 1
    assert bond3.src_sugar_index == 2
    assert bond3.dst_sugar_index == 3
    assert bond3.src_atom == 6
    assert bond3.dst_atom == 1
    assert bond4.src_sugar_index == 3
    assert bond4.dst_sugar_index == 4
    assert bond4.src_atom == 6
    assert bond4.dst_atom == 1


def test_complex_parsing_2():
    glycan = "MAN(4-1 FUC(4-1 MAN)(6-1 FUC(4-1 MAN)))(6-1 MAN(6-1 MAN(4-1 MAN)(6-1 FUC)))".replace(
        " ", ""
    )
    sugars, bonds = _glycan_string_to_sugars_and_bonds(glycan)
    assert len(sugars) == 9

    expected_bonds = [
        (0, 1),
        (1, 2),
        (1, 3),
        (3, 4),
        (0, 5),
        (5, 6),
        (6, 7),
        (6, 8),
    ]
    for (expected_src, expected_dst), bond in zip(expected_bonds, bonds, strict=True):
        assert bond.src_sugar_index == expected_src
        assert bond.dst_sugar_index == expected_dst


def test_glycan_tokenization_with_bond():
    """Test that tokenization works, and that atoms are dropped as expected."""
    glycan = ">glycan|foo\nNAG(4-1 NAG)\n"
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        fasta_file = tmp_path / "input.fasta"
        fasta_file.write_text(glycan)

        output_dir = tmp_path / "out"

        feature_context = make_all_atom_feature_context(
            fasta_file,
            output_dir=output_dir,
            use_esm_embeddings=False,  # Just a test; no need
        )

    # Each NAG component is C8 H15 N O6 -> 8 + 1 + 6 = 15 heavy atoms
    # The bond between them displaces one oxygen, leaving 2 * 15 - 1 = 29 atoms
    assert feature_context.structure_context.atom_exists_mask.sum() == 29
    # We originally constructed all atoms in dropped the atoms that leave
    assert feature_context.structure_context.atom_exists_mask.numel() == 30
    elements = Counter(
        feature_context.structure_context.atom_ref_element[
            feature_context.structure_context.atom_exists_mask
        ].tolist()
    )
    assert elements[6] == 16  # 6 = Carbon
    assert elements[7] == 2  # 7 = Nitrogen
    assert elements[8] == 11  # 8 = Oxygen

    # Single bond feature between O and C
    left, right = feature_context.structure_context.atom_covalent_bond_indices
    assert left.numel() == right.numel() == 1
    bond_elements = set(
        [
            feature_context.structure_context.atom_ref_element[left].item(),
            feature_context.structure_context.atom_ref_element[right].item(),
        ]
    )
    assert bond_elements == {8, 6}

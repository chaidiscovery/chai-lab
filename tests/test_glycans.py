# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
import pytest

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
    assert bond2.src_sugar_index == 0
    assert bond2.dst_sugar_index == 2
    assert bond3.src_sugar_index == 2
    assert bond3.src_sugar_index == 2
    assert bond3.dst_sugar_index == 3
    assert bond4.src_sugar_index == 3
    assert bond4.dst_sugar_index == 4


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

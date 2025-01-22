# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.


import pytest

from chai_lab.data.io.cif_utils import get_chain_letter


def test_get_chain_letter():
    with pytest.raises(AssertionError):
        get_chain_letter(0)
    assert get_chain_letter(1) == "A"
    assert get_chain_letter(26) == "Z"
    assert get_chain_letter(27) == "a"
    assert get_chain_letter(52) == "z"

    assert get_chain_letter(53) == "AA"
    assert get_chain_letter(54) == "AB"

    # For one-letter codes, there are 26 + 26 = 52 codes
    # For two-letter codes, there are 52 * 52 codes
    assert get_chain_letter(52 * 52 + 52) == "zz"

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
"""
Test for tokenization
"""

import numpy as np

from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order


def test_tokenization_basic():
    test_sequence = "RKDES"

    out, dels = tokenize_sequences_to_arrays([test_sequence])
    assert out.shape == dels.shape == (1, 5)
    assert np.all(
        out
        == np.array(
            [residue_types_with_nucleotides_order[res] for res in test_sequence]
        )
    )


def test_tokenization_with_insertion():
    """Insertions (lower case) should be ignored."""
    test_sequence = "RKDES"
    test_with_ins = "RKrkdesDES"

    out, dels = tokenize_sequences_to_arrays([test_sequence, test_with_ins])
    assert out.shape == dels.shape == (2, 5)
    assert np.all(out[0] == out[1])
    assert dels.sum() == 5
    assert dels[1, 2] == 5

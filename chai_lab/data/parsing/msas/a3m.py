# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Code for tokenizing aligned sequences as found in .a3m format.

Note that while we don't directly load from .a3m files for the model, the aligned.pqt format
we do use expects a column of sequences in the same format that .a3m gives.
"""

import re
import string
from functools import lru_cache
from io import StringIO
from pathlib import Path
from typing import Final

import numba
import numpy as np

from chai_lab.data.parsing.fasta import Fasta, read_fasta_content
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order

MAPPED_TOKEN_SKIP: Final[int] = -1
MAPPED_TOKEN_INSERTION: Final[int] = -2


@lru_cache
def _get_tokenization_mapping() -> np.ndarray:
    """
    each character in alignment either translated to token,
    or it should be ignored, or skipped but counted as insertion.

    tokenization_mapping allows doing this efficiently within numba,
     using just one lookup
    """
    # by default fill with MAPPED_TOKEN_SKIP
    mapping = np.full([256], dtype="int32", fill_value=MAPPED_TOKEN_SKIP)

    for letter in string.ascii_uppercase.encode():
        # init: all capital letters are treated as unknown ('X')
        # this will be changed for letters of aminoacids that we can process
        mapping[letter] = residue_types_with_nucleotides_order["X"]

    for letter_str, value in residue_types_with_nucleotides_order.items():
        if len(letter_str) == 1:  # ignore multi-letter
            mapping[ord(letter_str.encode())] = value

    for letter in string.ascii_lowercase.encode():
        mapping[letter] = MAPPED_TOKEN_INSERTION  # lower-case, insertion
    assert mapping[ord(b"-")] >= 0, "skip (-) should have its token"
    assert mapping[ord(b".")] == MAPPED_TOKEN_SKIP
    return mapping


@numba.jit(nopython=True)
def _parse_seqs_to_ndarrays(
    # \n added after every alignment; single alignment has no \n
    alignments_concatenated: bytes,
    tokenization_mapping: np.ndarray,
    # pre-allocated numpy arrays to be filled
    out_sequences: np.ndarray,
    out_deletions: np.ndarray,
):
    """Fills in out_sequences and out_deletions in-place."""
    start_strpos = 0
    seq_index = -1  # first used value is zero
    SPLITTER = 10  # ord(b"\n") = 10
    assert alignments_concatenated[-1] == SPLITTER
    n_positions = out_sequences.shape[1]
    for current_strpos, char in enumerate(alignments_concatenated):
        if char == SPLITTER:
            aligned_seq = alignments_concatenated[start_strpos:current_strpos]
            start_strpos = current_strpos + 1
            seq_index += 1

            # process single aligned sequence
            current_n_skips = 0
            position = 0
            for char_int in aligned_seq:
                mapped = tokenization_mapping[char_int]
                if mapped == MAPPED_TOKEN_SKIP:
                    pass
                elif mapped == MAPPED_TOKEN_INSERTION:
                    current_n_skips += 1
                else:
                    out_sequences[seq_index, position] = mapped
                    out_deletions[seq_index, position] = min(current_n_skips, 255)
                    current_n_skips = 0
                    position += 1

            # assert we don't need any padding; should be the case for valid sequences
            assert position == n_positions, (position, n_positions, seq_index)


def tokenize_sequences_to_arrays(
    seqs_str: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Tokenize a list of aligned sequences in a3m format."""
    # see: https://yanglab.qd.sdu.edu.cn/trRosetta/msa_format.html#a3m
    assert seqs_str, "Must provide non-empty list of sequences."
    seq_len: int = sum(c in string.ascii_uppercase or c == "-" for c in seqs_str[0])
    n_seqs = len(seqs_str)
    out_sequences = np.zeros(shape=(n_seqs, seq_len), dtype="uint8")
    out_deletions = np.zeros(shape=(n_seqs, seq_len), dtype="uint8")
    _parse_seqs_to_ndarrays(
        alignments_concatenated=b"".join(s.encode() + b"\n" for s in seqs_str),
        tokenization_mapping=_get_tokenization_mapping(),
        out_sequences=out_sequences,
        out_deletions=out_deletions,
    )
    return out_sequences, out_deletions


def read_colabfold_a3m(fname: Path) -> dict[str, list[Fasta]]:
    """Returns mapping of MSA hits per identifier in the given a3m file.

    The query line in each block of MSA hits is retained.
    """
    text = fname.read_text()
    retval: dict[str, list[Fasta]] = {}
    for block in text.split("\x00"):  # Splits on null byte
        if not block:
            continue
        strio = StringIO(block)
        hits = read_fasta_content(strio)
        assert len(hits) > 0

        query = hits[0].header
        assert re.match(r"^[0-9]{3}$", query)
        retval[query] = hits
    return retval

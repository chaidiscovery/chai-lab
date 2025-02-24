# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
from dataclasses import dataclass
from pathlib import Path

import torch
from torch import Tensor

from chai_lab.data import residue_constants as rc
from chai_lab.utils.typing import Bool, Int32, UInt8, typecheck


@typecheck
@dataclass(frozen=True)
class TemplateHit:
    """Template hit."""

    query_pdb_id: str
    query_sequence: str
    index: int
    pdb_id: str
    chain_id: str
    # The interval is half-open, i.e., [start, end)
    # These are basically counting the residues that aren't "-" gaps.
    hit_start: int  # Start is inclusive; within the full template sequece
    hit_end: int  # End is exclusive; within the full template sequence
    hit_tokens: Int32[Tensor, "n_tokens"]
    deletion_matrix: UInt8[Tensor, "n_tokens"]
    query_seq_realigned: str  # The realigned sequence
    cif_path: Path | None = None

    def __post_init__(self):
        q_i = self.indices_query
        h_i = self.indices_hit
        if q_i.shape != h_i.shape:
            raise ValueError(
                f"Query/hit indices mismatched size: {q_i.shape} {h_i.shape} - {q_i} {h_i}"
            )
        if h_i[0] != self.hit_start or h_i[-1] != self.hit_end - 1:
            raise ValueError(
                f"{self.query_pdb_id=} <> {self.pdb_id=} yields discrepant start/end "
                f"indices: {h_i[0]} {h_i[-1]} != {self.hit_start} {self.hit_end} - 1."
                "First set of indices is from parsing hit_tokens; second set is from "
                "hit metadata in source file."
            )

        if not self.hit_end > self.hit_start:
            raise ValueError(f"Invalid hit start/end: {self.hit_start} {self.hit_end}")

        if not self.indices_hit_within_subregion.min() >= 0:
            raise ValueError("negative indices")

        if self.cif_path is not None:
            assert self.cif_path.is_file()

    def __str__(self) -> str:
        fields = ", ".join(
            [
                f"query_pdb_id={self.query_pdb_id}",
                f"hit_pdb_id={self.pdb_id}",
                f"hit_chain={self.chain_id}",
                f"hit_start={self.hit_start}",
                f"hit_end={self.hit_end}",
                f"index={self.index}",
            ]
        )
        return f"TemplateHit({fields})"

    @property
    def hit_sequence(self) -> str:
        """Sequence of hit."""
        tokens = self.hit_tokens.tolist()
        return "".join([rc.residue_types_with_nucleotides[i] for i in tokens])

    @property
    def indices_query(self) -> Int32[Tensor, "m"]:
        """Indices of query residues corresponding to hit."""
        # Trim out the trailing gaps
        not_gap = self.hit_tokens != rc.residue_types_with_nucleotides_order["-"]
        matched_indices: list[int] = torch.where(not_gap)[0].tolist()
        # Returned indices are start inclusive end exclusive, i.e., [start, end)
        return torch.arange(start=matched_indices[0], end=matched_indices[-1] + 1)

    @property
    def indices_hit_within_subregion(self) -> Int32[Tensor, "n"]:
        """Indices of hit, accounting for ins/del in the hit alignment."""
        # NOTE These do NOT correspond 1:1 to the hit tokens! Because the hit tokens
        # have already been cropped down to match and these are the indices of the hit
        # BEFORE said cropping. These should however result in the hit_tokens when
        # applied to a fresh copy of the sequence.
        inc = torch.ones_like(self.hit_tokens, dtype=torch.int32)
        # Account for deletions that were dropped from the hit alignment
        # If there were N deletions to the left at i, then we need to forward the index
        # at i by N
        # If there is a gap in the hit, then we repeat the last index
        inc += self.deletion_matrix
        # Account for gaps in the hit alignment
        inc -= (self.hit_tokens == rc.residue_types_with_nucleotides_order["-"]).int()
        # Remove leading and trailing zero values
        nonzero_indices = (inc != 0).nonzero(as_tuple=True)[0]
        first_nonzero = nonzero_indices[0]
        last_nonzero = nonzero_indices[-1]
        retval = torch.cumsum(inc[first_nonzero : last_nonzero + 1], dim=0) - 1

        return retval[retval >= 0]

    @property
    def indices_hit(self) -> Int32[Tensor, "n"]:
        """Indices of hit within full hit sequence."""
        return self.indices_hit_within_subregion + self.hit_start

    @property
    def hit_valid_mask(self) -> Bool[Tensor, "n"]:
        """Mask of same size as indices_hit that is false for hit_indices positions that
        are gaps.

        Gaps here are like:
        query: RKDES    NOT     RKDES
        hit:   RK-ES    NOT     RKgDES
        This is because in the "-" gap instances, we repeat the last index in
        indices_hit to match on size, but we want to mask out these positions.
        """
        offset = self.indices_hit[1:] - self.indices_hit[:-1]
        # Gaps have 0 offset from previous; set these to False
        # Insertions have > 1 offset, but these are still "valid" so set to True
        retval = torch.cat((torch.tensor([True]), offset > 0))
        assert retval.shape == self.indices_hit.shape
        return retval

    @property
    def query_start_end(self) -> tuple[int, int]:
        start = int(self.indices_query.min().item())
        end = int(self.indices_query.max().item())
        return start, end

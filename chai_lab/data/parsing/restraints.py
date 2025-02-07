# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
"""
Parse table where each row correpsonds to a pairwise interaction
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
import pandera as pa
from pandera.typing import Series
from typing_extensions import assert_never

from chai_lab.utils.typing import typecheck


class PairwiseInteractionType(Enum):
    COVALENT = "covalent"
    CONTACT = "contact"
    POCKET = "pocket"


class PairwiseConstraintDataframeModel(pa.DataFrameModel):
    """Model for constraints table."""

    restraint_id: Series[str] = pa.Field(coerce=True, unique=True)
    chainA: Series[str] = pa.Field(nullable=False, coerce=True)
    res_idxA: Series[str] = pa.Field(nullable=True, coerce=True)
    chainB: Series[str] = pa.Field(nullable=False, coerce=True)
    res_idxB: Series[str] = pa.Field(nullable=True, coerce=True)
    max_distance_angstrom: Series[float] = pa.Field(nullable=True, ge=0.0)
    min_distance_angstrom: Series[float] = pa.Field(nullable=True, ge=0.0)
    connection_type: Series[str] = pa.Field(
        isin=[i.value for i in PairwiseInteractionType]
    )
    confidence: Series[float] = pa.Field(ge=0.0, le=1.0, nullable=True)
    comment: Series[str] = pa.Field(nullable=True)


@typecheck
@dataclass(frozen=True)
class PairwiseInteraction:
    """Interaction between two tokens. Indices are expected to be 1-based."""

    chainA: str
    res_idxA: str
    atom_nameA: str
    chainB: str
    res_idxB: str
    atom_nameB: str
    connection_type: PairwiseInteractionType
    max_dist_angstrom: float = 10.0
    min_dist_angstrom: float = 0.0
    confidence: float = 1.0
    comment: str = ""

    def __post_init__(self):
        # Chains should always be provided for both partners in an interaction
        assert (
            self.chainA and self.chainB
        ), f"Chains must be non-null: {self.chainA=} {self.chainB=}"
        if ~np.isnan(self.max_dist_angstrom) and ~np.isnan(self.min_dist_angstrom):
            assert self.max_dist_angstrom >= self.min_dist_angstrom

        match conn := self.connection_type:
            case PairwiseInteractionType.COVALENT:
                pass
            case PairwiseInteractionType.POCKET:
                assert self.res_idxA == "", "A (chain-level) should NOT specify a token"
                assert self.res_idxB != "", "B (token-level) should specify a token"
                assert (
                    self.atom_nameA == self.atom_nameB == ""
                ), "No atoms should be specified"
            case PairwiseInteractionType.CONTACT:
                assert self.res_idxA or self.atom_nameA, "A should specify a token/atom"
                assert self.res_idxB or self.atom_nameB, "B should specify a token/atom"
            case _:
                assert_never(conn)

        # Check that residue indices are > 0 (1-indexed)
        if self.res_idxA:
            assert self.res_idxA_pos > 0
        if self.res_idxB:
            assert self.res_idxB_pos > 0

    @property
    def res_idxA_name(self) -> str:
        """Single-char name of residue A."""
        return self.res_idxA[0] if self.res_idxA else ""

    @property
    def res_idxA_pos(self) -> int:
        """1-indexed position of residue A; defaults to 1 if not given."""
        # NOTE 1 default is because 1 is the minimum index under 1-indexing
        s = self.res_idxA[1:]
        return int(s) if s else 1

    @property
    def res_idxB_name(self) -> str:
        """Single-char name of residue B."""
        return self.res_idxB[0] if self.res_idxB else ""

    @property
    def res_idxB_pos(self) -> int:
        """1-indexed position of residue B; defaults to 1 if not given."""
        s = self.res_idxB[1:]
        return int(s) if s else 1

    def to_table_entry(self) -> dict[str, str | float]:
        """Format as table entry, sans leading restraint_id column."""
        values = {
            "chainA": self.chainA,
            "res_idxA": (
                self.res_idxA + (("@" + self.atom_nameA) if self.atom_nameA else "")
            ),
            "chainB": self.chainB,
            "res_idxB": (
                self.res_idxB + (("@" + self.atom_nameB) if self.atom_nameB else "")
            ),
            "connection_type": self.connection_type.value,
            "confidence": self.confidence,
            "min_distance_angstrom": self.min_dist_angstrom,
            "max_distance_angstrom": self.max_dist_angstrom,
            "comment": self.comment,
        }
        return values


def _parse_res_idx(res_idx: str) -> tuple[str, str]:
    """Parse a residue index string.

    >>> _parse_res_idx("A219")
    ('A219', '')
    >>> _parse_res_idx("A219@CA")
    ('A219', 'CA')
    >>> _parse_res_idx("@C2")
    ('', 'C2')
    """
    if res_idx.endswith("@"):
        raise ValueError(f"Invalid residue index: {res_idx}")
    parts = res_idx.split("@")
    if len(parts) == 1:
        parts.append("")
    elif len(parts) > 2 or all(len(part) == 0 for part in parts):
        raise ValueError(f"Invalid residue index: {res_idx}")
    res, atomname = parts
    return res, atomname


def _parse_row(row: pd.Series) -> PairwiseInteraction:
    """Parse a row of a pairwise interaction table."""
    resA, atomA = _parse_res_idx(row["res_idxA"])
    resB, atomB = _parse_res_idx(row["res_idxB"])
    return PairwiseInteraction(
        chainA=row["chainA"],
        res_idxA=resA,
        atom_nameA=atomA,
        chainB=row["chainB"],
        res_idxB=resB,
        atom_nameB=atomB,
        max_dist_angstrom=row["max_distance_angstrom"],
        min_dist_angstrom=row["min_distance_angstrom"],
        connection_type=PairwiseInteractionType(row["connection_type"]),
        confidence=row["confidence"],
        comment=row["comment"],
    )


def parse_pairwise_table(table: str | Path) -> list[PairwiseInteraction]:
    """Parse a table of pairwise interactions."""
    df = pd.read_csv(table)

    PairwiseConstraintDataframeModel.validate(df, inplace=True)

    # We only have floats and strings; this applies to strings
    df.update(df.select_dtypes(include=["object"]).fillna(""))

    # Set confidence to 1.0 if not provided
    df["confidence"] = df["confidence"].fillna(1.0)

    # Parsing
    return [_parse_row(row) for _, row in df.iterrows()]


def write_pairwise_table(interactions: list[PairwiseInteraction], fname: str | Path):
    """Write the interactions to a .csv file."""
    entries = []
    for i, entry in enumerate([i.to_table_entry() for i in interactions]):
        entry["restraint_id"] = f"restraint_{i}"
        entries.append(entry)
    df = pd.DataFrame.from_records(entries)

    df.to_csv(fname, index=False)
    return


if __name__ == "__main__":
    import doctest

    doctest.testmod()

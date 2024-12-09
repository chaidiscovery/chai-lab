# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Parsing logic for glycans.

- Each sugar in the glycan gets a distinct residue
    - Each residue has an increasing label seq

Example glycan strings:
MAN(6-1 FUC)(4-1 MAN)
"""

import re
from functools import lru_cache

from attr import dataclass

from chai_lab.data import residue_constants as rc
from chai_lab.data.parsing.structure.residue import Residue


@dataclass(frozen=True)
class GlycosidicBond:
    src_sugar_index: int  # 0-indexed
    dst_sugar_index: int  # 0-indexed
    src_atom: int  # 1-indexed
    dst_atom: int  # 1-indexed

    def __post_init__(self):
        assert self.src_sugar_index != self.dst_sugar_index
        assert self.src_atom > 0 and self.dst_atom > 0

    @property
    def src_atom_name(self) -> str:
        """Links between sugars are O-glycosidic bonds; we use src O dst C."""
        return f"O{self.src_atom}"

    @property
    def dst_atom_name(self) -> str:
        """Links between sugars are O-glycosidic bonds; we use src O dst C."""
        return f"C{self.dst_atom}"


@lru_cache(maxsize=32)
def _glycan_string_to_sugars_and_bonds(
    glycan_string: str,
) -> tuple[list[str], list[GlycosidicBond]]:
    """Parses the glycan string to its constituent sugars and bonds."""
    glycan_string = glycan_string.strip()  # Remove leading/trailing spaces
    sugars: list[str] = []  # Tracks all sugars
    parent_sugar_idx: list[int] = []  # Tracks the parent sugar for bond formation
    bonds: list[GlycosidicBond] = []
    open_count, closed_count = 0, 0
    i = 0  # We increment unevenly so manually handle
    while i < len(glycan_string):
        char = glycan_string[i]
        if char == " ":  # Space; skip
            i += 1
            continue
        if char == "(":  # Open bracket
            i += 1
            open_count += 1
            continue
        if char == ")":  # Close bracket
            closed_count += 1
            parent_sugar_idx.pop()  # Remove
            i += 1
            continue
        # Not a bracket or a space - should be either bond info or CCD
        chunk = glycan_string[i : i + 3]
        if re.match(r"[1-6]{1}-[1-6]{1}", chunk):
            s, d = chunk.split("-")
            assert parent_sugar_idx
            bonds.append(
                GlycosidicBond(
                    src_sugar_index=parent_sugar_idx[-1],
                    dst_sugar_index=len(sugars),  # Anticipate next
                    src_atom=int(s),
                    dst_atom=int(d),
                )
            )
            i += 3
        elif re.match(r"[0-9A-Z]{3}", chunk):
            sugars.append(chunk)
            parent_sugar_idx.append(len(sugars) - 1)  # latest sugar
            i += 3
        else:
            raise ValueError(f"Invalid glycan string: {glycan_string}")
    assert open_count == closed_count
    return sugars, bonds


def glycan_string_residues(glycan_string: str) -> list[Residue]:
    sugars, _bonds = _glycan_string_to_sugars_and_bonds(glycan_string)
    if not sugars:
        raise ValueError(f"No residues parsed from {glycan_string=}")
    return [
        Residue(
            name=sugar,
            label_seq=i + 1,
            restype=rc.residue_types_with_nucleotides_order["X"],
            residue_index=i,
            is_missing=False,
            b_factor_or_plddt=0.0,
            conformer_data=None,
            is_covalent_bonded=True,
        )
        for i, sugar in enumerate(sugars)
    ]

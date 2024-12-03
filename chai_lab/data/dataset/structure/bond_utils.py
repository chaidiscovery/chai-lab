# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from torch import Tensor

from chai_lab.data.parsing.glycans import _glycan_string_to_sugars_and_bonds
from chai_lab.utils.typing import Int, typecheck


@typecheck
def get_atom_covalent_bond_pairs_from_glycan_string(
    glycan_string: str,
    token_residue_index: Int[Tensor, "n_tokens"],
    atom_ref_name: list[str],
) -> tuple[Int[Tensor, "n_bonds"], Int[Tensor, "n_bonds"]]:
    """Infer bond pairs between glycans sugar rings."""
    if not glycan_string:
        return (
            torch.zeros(0, dtype=torch.long),
            torch.zeros(0, dtype=torch.long),
        )

    assert token_residue_index.numel() == len(atom_ref_name)
    _sugars, bonds = _glycan_string_to_sugars_and_bonds(glycan_string)
    left_bonds, right_bonds = [], []
    for bond in bonds:
        left_chain_mask = token_residue_index == bond.src_sugar_index
        left_res_mask = [n == bond.src_atom_name for n in atom_ref_name]
        right_chain_mask = token_residue_index == bond.dst_sugar_index
        right_res_mask = [n == bond.dst_atom_name for n in atom_ref_name]
        left_res = left_chain_mask & torch.tensor(left_res_mask)
        right_res = right_chain_mask & torch.tensor(right_res_mask)
        assert left_res.sum() == 1, f"Expected unique atom, got {left_res=}"
        assert right_res.sum() == 1, f"Expected unique atom, got {right_res=}"
        left_res_idx, *_ = torch.where(left_res)
        right_res_idx, *_ = torch.where(right_res)
        left_bonds.append(left_res_idx.item())
        right_bonds.append(right_res_idx.item())
    bonds_to_add = (
        torch.tensor(left_bonds, dtype=torch.int),
        torch.tensor(right_bonds, dtype=torch.int),
    )
    return bonds_to_add

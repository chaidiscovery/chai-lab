# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.


import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data import residue_constants as rc
from chai_lab.data.parsing.glycans import _glycan_string_to_sugars_and_bonds
from chai_lab.data.parsing.restraints import (
    PairwiseInteraction,
    PairwiseInteractionType,
)
from chai_lab.model.utils import get_asym_id_from_subchain_id
from chai_lab.utils.tensor_utils import string_to_tensorcode
from chai_lab.utils.typing import Int, UInt8, typecheck


@typecheck
def get_atom_covalent_bond_pairs_from_constraints(
    provided_constraints: list[PairwiseInteraction],
    token_residue_index: Int[Tensor, "n_tokens"],
    token_residue_name: UInt8[Tensor, "n_tokens 8"],
    token_subchain_id: UInt8[Tensor, "n_tokens 4"],
    token_asym_id: Int[Tensor, "n_tokens"],
    atom_token_index: Int[Tensor, "n_atoms"],
    atom_ref_name: list[str],
) -> tuple[Int[Tensor, "n_bonds"], Int[Tensor, "n_bonds"]]:
    """Determine bond pairs, which are returned as atom indices that are bonded."""
    ret_a: list[int] = []
    ret_b: list[int] = []
    for constraint in provided_constraints:
        match ctype := constraint.connection_type:
            case PairwiseInteractionType.COVALENT:
                assert (
                    constraint.atom_nameA and constraint.atom_nameB
                ), "Atoms must be provided for covalent bonds"
                # Figure out the asym id that we care about
                left_asym_id = get_asym_id_from_subchain_id(
                    subchain_id=constraint.chainA,
                    source_pdb_chain_id=token_subchain_id,
                    token_asym_id=token_asym_id,
                )
                left_token_asym_mask = token_asym_id == left_asym_id
                right_asym_id = get_asym_id_from_subchain_id(
                    subchain_id=constraint.chainB,
                    source_pdb_chain_id=token_subchain_id,
                    token_asym_id=token_asym_id,
                )
                right_token_asym_mask = token_asym_id == right_asym_id
                assert torch.any(left_token_asym_mask) and torch.any(
                    right_token_asym_mask
                )

                # Get the token index that we want
                left_token_index_mask = (
                    token_residue_index == constraint.res_idxA_pos - 1
                )
                right_token_index_mask = (
                    token_residue_index == constraint.res_idxB_pos - 1
                )
                assert torch.any(left_token_index_mask) and torch.any(
                    right_token_index_mask
                )

                # Combine these to get the specific residue specified
                left_residue_mask = left_token_asym_mask & left_token_index_mask
                if constraint.res_idxA_name:
                    three_letter = string_to_tensorcode(
                        rc.restype_1to3.get(constraint.res_idxA_name, "UNK"),
                        pad_to_length=token_residue_name.shape[-1],
                    )
                    resname_matches = (
                        token_residue_name == rearrange(three_letter, "d -> 1 d")
                    ).all(dim=-1)
                    assert resname_matches.shape == left_residue_mask.shape
                    left_residue_mask &= resname_matches
                right_residue_mask = right_token_asym_mask & right_token_index_mask
                if constraint.res_idxB_name:
                    three_letter = string_to_tensorcode(
                        rc.restype_1to3.get(constraint.res_idxB_name, "UNK"),
                        pad_to_length=token_residue_name.shape[-1],
                    )
                    resname_matches = (
                        token_residue_name == rearrange(three_letter, "d -> 1 d")
                    ).all(dim=-1)
                    assert resname_matches.shape == right_residue_mask.shape
                    right_residue_mask &= resname_matches
                # NOTE there are multiple residues in these residue masks due to
                # per-atom tokenization of glycans
                # These indices do not reset for new chains (matching atom_token_index)
                left_residue_idx = torch.where(left_residue_mask)[0]
                right_residue_idx = torch.where(right_residue_mask)[0]
                assert left_residue_idx.numel() > 0 and right_residue_idx.numel() > 0

                # Find the atoms belonging to this residue
                left_atoms_mask = torch.isin(
                    atom_token_index, test_elements=left_residue_idx
                )
                right_atoms_mask = torch.isin(
                    atom_token_index, test_elements=right_residue_idx
                )
                assert torch.any(left_atoms_mask) and torch.any(right_atoms_mask)

                # Find atoms matching on atom name
                left_name_mask = torch.tensor(
                    [n == constraint.atom_nameA for n in atom_ref_name],
                    dtype=torch.bool,
                )
                right_name_mask = torch.tensor(
                    [n == constraint.atom_nameB for n in atom_ref_name],
                    dtype=torch.bool,
                )

                left_atom_mask = left_atoms_mask & left_name_mask
                right_atom_mask = right_atoms_mask & right_name_mask
                assert (
                    torch.sum(left_atom_mask) == torch.sum(right_atom_mask) == 1
                ), f"Expect single atoms, got {torch.sum(left_atom_mask)}, {torch.sum(right_atom_mask)}"

                (left_atom_idx,) = torch.where(left_atom_mask)
                (right_atom_idx,) = torch.where(right_atom_mask)
                ret_a.append(left_atom_idx.item())  # type: ignore
                ret_b.append(right_atom_idx.item())  # type: ignore

            case PairwiseInteractionType.CONTACT | PairwiseInteractionType.POCKET:
                # These are handled as constraints, not as bonds
                pass
            case _:
                raise ValueError(f"Unrecognized pariwise interaction: {ctype}")
    return torch.tensor(ret_a, dtype=torch.long), torch.tensor(ret_b, dtype=torch.long)


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

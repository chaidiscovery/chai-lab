# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.data.features.token_utils import get_centre_positions_and_mask
from chai_lab.utils.tensor_utils import cdist, und_self
from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
def abc_is_colinear(
    atoms_a: Float[Tensor, "b n_triplets 3"],
    atoms_b: Float[Tensor, "b n_triplets 3"],
    atoms_c: Float[Tensor, "b n_triplets 3"],
) -> Bool[Tensor, "b n_triplets"]:
    """Check to see if each triplet of 3 atoms (a, b, c) is co-linear."""
    w1 = atoms_a - atoms_b
    w1 /= torch.linalg.norm(w1, dim=-1, keepdim=True)
    w2 = atoms_c - atoms_b
    w2 /= torch.linalg.norm(w2, dim=-1, keepdim=True)

    cos_sim = torch.sum(w1 * w2, dim=-1)
    cos_sim = torch.clamp(cos_sim, -1.0, 1.0)
    angle = torch.acos(cos_sim)  # radians

    # Colinearity should cover cases that are very small acute angles and cases of large
    # obtuse angles that are close to 180 degrees.
    colinear = (
        torch.isnan(angle)
        | (angle < 25 / 180 * torch.pi)
        | (angle > 155 / 180 * torch.pi)
    )
    return colinear


@typecheck
def get_single_atom_frames(
    atom_coords: Float[Tensor, "b n_atoms 3"],
    token_asym_id: Int[Tensor, "b n_tokens"],
    token_residue_index: Int[Tensor, "b n_tokens"],
    token_backbone_frame_mask: Bool[Tensor, "b n_tokens"],
    token_centre_atom_index: Int[Tensor, "b n_tokens"],
    token_exists_mask: Bool[Tensor, "b n_tokens"],
    atom_exists_mask: Bool[Tensor, "b n_atoms"],
    atom_token_index: Int[Tensor, "b n_atoms"],
) -> tuple[Int[Tensor, "b n_tokens 3"], Bool[Tensor, "b n_tokens"]]:
    """Makes frames for everything that is tokenized per-atom"""
    # For tokens that are one atom per token, a_i, b_i, c_i for frame is:
    # - token atom is assigned as b_i
    # - closest atom to the token atom is a_i
    # - second closest atom is c_i

    # Compute distances; n_tokens size
    centre_coords, centre_mask = get_centre_positions_and_mask(
        atom_coords,
        atom_exists_mask,
        token_centre_atom_index,
        token_exists_mask,
    )

    asym_mask = rearrange(token_asym_id, "b i -> b i 1") == rearrange(
        token_asym_id, "b j -> b 1 j"
    )
    res_idx_mask = rearrange(token_residue_index, "b i -> b i 1") == rearrange(
        token_residue_index, "b j -> b 1 j"
    )
    dists = cdist(centre_coords)  # Symmetric (tokens x tokens)
    # Mask out distances that don't exist
    centre_mask_square = und_self(centre_mask, "b i, b j -> b i j")
    # restrict to intra-residue pairs with valid coords
    dists = dists.masked_fill(
        ~centre_mask_square | ~asym_mask | ~res_idx_mask, torch.inf
    )

    B, tokens = dists.shape[:2]
    device = dists.device

    # Mask out diagonal
    batch_indices = torch.arange(B, device=device)[..., None, None]
    dists[batch_indices, torch.eye(tokens, device=device).bool()] = torch.inf

    _, idces = torch.topk(dists, 2, dim=-1, largest=False)  # b, n_tokens, 2
    a, c = idces.unbind(dim=-1)
    b = torch.arange(tokens, device=device).unsqueeze(0)  # Token index

    # Convert from token index to ATOM index
    batch_indices = torch.arange(B, device=device)[..., None]
    abc_atom_indices = torch.stack(
        [token_centre_atom_index[batch_indices, idx] for idx in [a, b, c]],
        dim=-1,
    )
    abc_coords_mask = torch.stack(
        [centre_mask[batch_indices, idx] for idx in [a, b, c]],
        dim=-1,
    ).all(dim=-1)

    # Make mask for tokens within the same chain
    a_res_idx = token_residue_index[batch_indices, a]
    b_res_idx = token_residue_index[batch_indices, b]
    c_res_idx = token_residue_index[batch_indices, c]

    a_asym, b_asym, c_asym = (
        token_asym_id[batch_indices, a],
        token_asym_id[batch_indices, b],
        token_asym_id[batch_indices, c],
    )

    same_residue = (a_res_idx == b_res_idx) & (b_res_idx == c_res_idx)
    same_chain = (a_asym == b_asym) & (b_asym == c_asym)

    # Check for co-linearity (< 25 degrees deviation)
    colinear = abc_is_colinear(
        centre_coords[batch_indices, a],
        centre_coords[batch_indices, b],
        centre_coords[batch_indices, c],
    )

    # Positions where the token backbone was NOT already defined, shares the same
    # entity_id, are not co-linear, and is actually a centre atom
    mask = torch.ones_like(token_backbone_frame_mask)
    for i in range(mask.shape[0]):
        all_idces, counts = torch.unique(atom_token_index[i], return_counts=True)
        not_single_idces = all_idces[counts != 1]
        mask[i, not_single_idces] = False

    mask &= (
        ~token_backbone_frame_mask
        & same_residue
        & same_chain
        & ~colinear
        & abc_coords_mask
        & token_exists_mask
    )

    return abc_atom_indices, mask


@typecheck
def get_frames_and_mask(
    atom_coords: Float[Tensor, "b n_atoms 3"],
    token_asym_id: Int[Tensor, "b n_tokens"],
    token_residue_index: Int[Tensor, "b n_tokens"],
    token_backbone_frame_mask: Bool[Tensor, "b n_tokens"],
    token_centre_atom_index: Int[Tensor, "b n_tokens"],
    token_exists_mask: Bool[Tensor, "b n_tokens"],
    atom_exists_mask: Bool[Tensor, "b n_atoms"],
    backbone_frame_idces: Int[Tensor, "b n_tokens 3"],
    atom_token_index: Int[Tensor, "b n_atoms"],
) -> tuple[Int[Tensor, "b n_tokens 3"], Bool[Tensor, "b n_tokens"]]:
    """Computes union of defined backbone frames and single atom frames"""
    single_atom_frame_idces, single_atom_frames_mask = get_single_atom_frames(
        atom_coords=atom_coords,
        token_asym_id=token_asym_id,
        token_residue_index=token_residue_index,
        token_backbone_frame_mask=token_backbone_frame_mask,
        token_centre_atom_index=token_centre_atom_index,
        token_exists_mask=token_exists_mask,
        atom_exists_mask=atom_exists_mask,
        atom_token_index=atom_token_index,
    )

    frame_idces = backbone_frame_idces.clone()
    mask = repeat(single_atom_frames_mask, "b n -> b n 3")
    frame_idces[mask] = single_atom_frame_idces[mask]

    all_frames_mask = single_atom_frames_mask | token_backbone_frame_mask

    return frame_idces, all_frames_mask

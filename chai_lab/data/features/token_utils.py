# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import repeat
from torch import Tensor

from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
def get_centre_positions_and_mask(
    atom_gt_coords: Float[Tensor, "... n_atoms 3"],
    atom_exists_mask: Bool[Tensor, "... n_atoms"],
    token_centre_atom_index: Int[Tensor, "... n_tokens"],
    token_exists_mask: Bool[Tensor, "... n_tokens"],
) -> tuple[Float[Tensor, "... n_tokens 3"], Bool[Tensor, "... n_tokens"]]:
    assert token_centre_atom_index.dtype in (torch.int32, torch.long)
    center_index = token_centre_atom_index.long()
    indices = repeat(center_index, "... n -> ... n c", c=3)
    center_pos = torch.gather(atom_gt_coords, dim=-2, index=indices)
    center_mask = torch.gather(atom_exists_mask, dim=-1, index=center_index)

    # because token_centre_atom_index is zero-padded, and because
    # atom number 0 is probably a valid atom, we need to reapply
    # the token mask
    center_mask = center_mask & token_exists_mask

    return center_pos, center_mask


@typecheck
def get_token_reference_atom_positions_and_mask(
    atom_pos: Float[Tensor, "b a d"],
    atom_mask: Bool[Tensor, "b a"],
    token_reference_atom_index: Int[Tensor, "b n"],
    token_exists_mask: Bool[Tensor, "b n_tokens"],
) -> tuple[Float[Tensor, "b n d"], Bool[Tensor, "b n"]]:
    """
    Get positions of the reference atom (CB for proteins, C2 or C4 for nucleic acids, sole atom for atom tokens)
    for each token
    """
    batch_indices = torch.arange(0, atom_pos.shape[0])[:, None]

    reference_atom_pos = atom_pos[batch_indices, token_reference_atom_index]
    reference_atom_mask = atom_mask[batch_indices, token_reference_atom_index]

    # because token_reference_atom_index is zero-padded, and because
    # atom number 0 is probably a valid atom, we need to reapply
    # the token mask
    reference_atom_mask = reference_atom_mask & token_exists_mask

    return reference_atom_pos, reference_atom_mask

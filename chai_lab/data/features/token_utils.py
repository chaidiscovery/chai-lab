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

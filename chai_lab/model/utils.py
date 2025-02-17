# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Any

import torch
from einops import rearrange, reduce, repeat
from torch import Tensor

from chai_lab.utils.tensor_utils import string_to_tensorcode, und
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck


def get_qkv_indices_for_blocks(
    sequence_length: int,
    stride: int,
    kv_block_size: int,
    device: Any,
) -> tuple[
    Int[torch.Tensor, "bl bl_q"],
    Int[torch.Tensor, "bl bl_kv"],
    Bool[torch.Tensor, "bl bl_kv"],
]:
    """Gets q, kv indices for local attention blocks."""
    sequence_length
    # from now on pretend q and kv are different axes
    num_blocks = sequence_length // stride
    assert (
        sequence_length == num_blocks * stride
    ), f"only seqlens divisible by {stride=} are supported, not {sequence_length=}"
    q_indices = torch.arange(sequence_length, device=device)
    q_indices = rearrange(
        q_indices, "(bl bl_q) -> bl bl_q", bl=num_blocks, bl_q=stride
    )  # bl bl_q -> q
    kv_indices = q_indices[:, :1] + (stride - kv_block_size) // 2
    kv_indices = kv_indices + torch.arange(
        kv_block_size, device=kv_indices.device
    )  # bl bl_kv -> kv
    # mask out positions where kv_indices gets wrapped
    # Rationale: the local attention block should always process
    # local blocks (i.e. same rel-positional encodings for each block.)
    kv_mask = (kv_indices < sequence_length) & (kv_indices >= 0)
    # Use of % not .clamp is important for short sequences
    kv_indices = kv_indices % sequence_length
    # q_idx is returned for reference, downstream code uses reshapes instead
    return q_indices, kv_indices, kv_mask


@typecheck
def get_block_atom_pair_mask(
    atom_single_mask: Bool[Tensor, "b a"],
    q_idx: Int[Tensor, "bl bl_q"],
    kv_idx: Int[Tensor, "bl bl_kv"],
    kv_is_wrapped_mask: Bool[Tensor, "bl bl_kv"],
) -> Bool[Tensor, "b bl bl_q bl_kv"]:
    atom_q_mask = atom_single_mask[:, q_idx]
    atom_kv_mask = atom_single_mask[:, kv_idx]

    block_atom_pair_mask = und(
        atom_q_mask, atom_kv_mask, "b bl bl_q, b bl bl_kv -> b bl bl_q bl_kv"
    )

    block_atom_pair_mask &= rearrange(kv_is_wrapped_mask, "bl bl_kv -> 1 bl 1 bl_kv")
    return block_atom_pair_mask


@typecheck
def calc_centroid(
    coords: Float[Tensor, "b a 3"],
    mask: Bool[Tensor, "#b a"],
    weights: Float[Tensor, "b a"] | None = None,
) -> Float[Tensor, "b 3"]:
    # mean-center coordinates
    masked_weights = weights * mask if weights is not None else mask.to(coords.dtype)
    masked_weights /= reduce(masked_weights, "b a -> b 1", "sum").clamp(min=1e-4)
    # not using einsum to avoid autocasting
    return reduce(coords * masked_weights[:, :, None], "b a c -> b c", "sum")


def _copysign(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """
    Transform from: https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms
    Return a tensor where each element has the absolute value taken from the,
    corresponding element of a, with sign taken from the corresponding
    element of b. This is like the standard copysign floating-point operation,
    but is not careful about negative 0 and NaN.

    Args:
        a: source tensor.
        b: tensor whose signs will be used, of the same shape as a.

    Returns:
        Tensor of the same shape as a with the signs of b.
    """
    signs_differ = (a < 0) != (b < 0)
    return torch.where(signs_differ, -a, a)


def quaternion_to_matrix(quaternions: torch.Tensor) -> torch.Tensor:
    """
    Transform from: https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms
    Convert rotations given as quaternions to rotation matrices.

    Args:
        quaternions: quaternions with real part first,
            as tensor of shape (..., 4).

    Returns:
        Rotation matrices as tensor of shape (..., 3, 3).
    """
    r, i, j, k = torch.unbind(quaternions, -1)
    # pyre-fixme[58]: `/` is not supported for operand types `float` and `Tensor`.
    two_s = 2.0 / (quaternions * quaternions).sum(-1)

    o = torch.stack(
        (
            1 - two_s * (j * j + k * k),
            two_s * (i * j - k * r),
            two_s * (i * k + j * r),
            two_s * (i * j + k * r),
            1 - two_s * (i * i + k * k),
            two_s * (j * k - i * r),
            two_s * (i * k - j * r),
            two_s * (j * k + i * r),
            1 - two_s * (i * i + j * j),
        ),
        -1,
    )
    return o.reshape(quaternions.shape[:-1] + (3, 3))


def random_quaternions(
    n: int, dtype: torch.dtype | None = None, device: Any | None = None
) -> torch.Tensor:
    """
    Transform from: https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms
    Generate random quaternions representing rotations,
    i.e. versors with nonnegative real part.

    Args:
        n: Number of quaternions in a batch to return.
        dtype: Type to return.
        device: Desired device of returned tensor. Default:
            uses the current device for the default tensor type.

    Returns:
        Quaternions as tensor of shape (N, 4).
    """
    if isinstance(device, str):
        device = torch.device(device)
    o = torch.randn((n, 4), dtype=dtype, device=device)
    s = (o * o).sum(1)
    o = o / _copysign(torch.sqrt(s), o[:, 0])[:, None]
    return o


def random_rotations(
    n: int, dtype: torch.dtype | None = None, device: Any = None
) -> torch.Tensor:
    """
    Transform from: https://pytorch3d.readthedocs.io/en/latest/_modules/pytorch3d/transforms
    Generate random rotations as 3x3 rotation matrices.

    Args:
        n: Number of rotation matrices in a batch to return.
        dtype: Type to return.
        device: Device of returned tensor. Default: if None,
            uses the current device for the default tensor type.

    Returns:
        Rotation matrices as tensor of shape (n, 3, 3).
    """
    quaternions = random_quaternions(n, dtype=dtype, device=device)
    return quaternion_to_matrix(quaternions)


@torch.no_grad()
@typecheck
def center_random_augmentation(
    atom_coords: Float[Tensor, "b a 3"],
    atom_single_mask: Bool[Tensor, "#b a"],
    s_trans: float = 1.0,
    rotations: Float[Tensor, "b 3 3"] | None = None,
) -> Float[Tensor, "b a 3"]:
    centroid = calc_centroid(atom_coords, mask=atom_single_mask)
    centroid = rearrange(centroid, "b c -> b 1 c")
    atom_coords = atom_coords - centroid
    # randomly rotate
    if rotations is None:
        rotations = random_rotations(atom_coords.shape[0], device=atom_coords.device)
    rotated_coords = torch.einsum("b i j, b a j -> b a i", rotations, atom_coords)
    random_translation = torch.randn_like(centroid)  # b 1 c=3
    return rotated_coords + s_trans * random_translation


@typecheck
def get_asym_id_from_subchain_id(
    subchain_id: str,
    source_pdb_chain_id: UInt8[Tensor, "n_tokens 4"],
    token_asym_id: Int[Tensor, "n"],
) -> int:
    # encde the subchain ids and perform lookup in context features
    chain_id_tensorcode = string_to_tensorcode(subchain_id, pad_to_length=4)
    chain_id_tensorcode = chain_id_tensorcode.to(token_asym_id.device)
    # create masks
    chain_id_tensorcode = repeat(chain_id_tensorcode, "c -> 1 c")
    chain_id_mask = torch.all(chain_id_tensorcode == source_pdb_chain_id, dim=-1)
    # check uniqueness
    chain_id_asyms = torch.unique(token_asym_id[chain_id_mask])

    assert len(chain_id_asyms) == 1, (
        f"Expected only one token asym, but got {len(chain_id_asyms)} "
        f"asyms: {chain_id_asyms}"
    )
    return chain_id_asyms[0].item()

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import typing
from dataclasses import asdict, is_dataclass, replace
from functools import lru_cache
from typing import TypeVar

import torch
import torch.nn.functional as F
from einops import rearrange
from torch import Tensor

from chai_lab.utils.defaults import default
from chai_lab.utils.typing import Bool, Float, UInt8, typecheck


@typecheck
def cdist(
    x: Float[Tensor, "... p m"],
    y: Float[Tensor, "... r m"] | None = None,
    p: float = 2.0,
) -> Float[Tensor, "... p r"]:
    y = default(y, x)
    assert x.ndim == y.ndim

    _threshold = 2147400000
    n, m = x.shape[-2], y.shape[-2]

    flat_size = torch.prod(torch.tensor(x.shape[:-2])) * n * m

    if x.is_cuda and flat_size > _threshold:
        # Torch cdist without mm fails when the total number of
        # points is > _threshold (in dimension 3)
        # or 8192 points for batch size 32.
        # To preserve accuracy, we fallback to naive distances
        return _naive_pairwise_distances(x, y)

    return torch.cdist(x1=x, x2=y, compute_mode="donot_use_mm_for_euclid_dist", p=p)


@typecheck
def _naive_pairwise_distances(
    x: Float[Tensor, "... p m"],
    y: Float[Tensor, "... r m"] | None = None,
    eps: float = 1e-10,
) -> Float[Tensor, "... p r"]:
    y = default(y, x)
    diff = x.unsqueeze(-2) - y.unsqueeze(-3)

    return diff.pow_(2).sum(dim=-1).add_(eps).sqrt_()


@typecheck
def masked_mean(
    mask: Bool[Tensor, "..."],
    value: Tensor,
    dim: int | tuple[int, ...],
    keepdim=False,
) -> Tensor:
    mask = mask.expand(*value.shape)
    num = torch.sum(mask * value, dim=dim, keepdim=keepdim)
    denom = torch.sum(mask, dim=dim, keepdim=keepdim).clamp(min=1)
    return num / denom


@typecheck
def one_hot(x: Tensor, v_bins: Tensor) -> Tensor:
    """One hot encoding; v_bins should N-1 bins where N is desired bins."""
    bins = torch.searchsorted(v_bins, x)
    return F.one_hot(bins, v_bins.shape[-1] + 1).float()


@lru_cache()
def _get_individual_und_patterns(multipattern: str) -> list[str]:
    assert isinstance(multipattern, str), "pattern goes as last argument"
    left_parts, right_part = multipattern.split("->")
    assert "(" not in right_part, "parenthesis not supported for now"
    result = []

    all_left_ids = set()
    all_left_parts_have_ellipsis = True

    for left_part in left_parts.split(","):
        left_ids = set(left_part.split())
        if "..." not in left_ids:
            all_left_parts_have_ellipsis = False
        all_left_ids.update(left_ids)
        right_parts = []
        for token in right_part.split():
            if token == "1" or token in left_ids:  # '...' should be in left ids
                right_parts.append(token)
            elif token.isidentifier():
                right_parts.append("1")
            elif token == "...":
                raise RuntimeError(
                    f"Ellipis not in one of left sides of {multipattern=}"
                )
            else:
                raise RuntimeError(f"Unknown {token=} in {multipattern=}")
        result.append(f"{left_part} -> " + " ".join(right_parts))

    if "..." in right_part.split():
        msg = "for now ALL or NONE left parts should have ellipsis (...) "
        assert all_left_parts_have_ellipsis, msg

    unk_ids = [
        x
        for x in right_part.split()
        if x not in all_left_ids and x != "1" and x != "..."
    ]
    assert len(unk_ids) == 0, f"{unk_ids=} not found on left side of {multipattern}"

    return result


@typing.overload
def und(t1: Tensor, pattern: str) -> Tensor: ...


@typing.overload
def und(t1: Tensor, t2: Tensor, pattern: str) -> Tensor: ...


@typing.overload
def und(t1: Tensor, t2: Tensor, t3: Tensor, pattern: str) -> Tensor: ...


@typing.overload
def und(t1: Tensor, t2: Tensor, t3: Tensor, t4: Tensor, pattern: str) -> Tensor: ...


def und(*args):
    """
    Micro-extension to einops.

    Performs & (logical_and) for several masks.
    Similar to einsum over masks, but additionally can add/remove 1-dims.

    > und(mask1, mask2, "b i, b j -> b 1 i j")
    """
    *tensors, multipattern = args
    patterns = _get_individual_und_patterns(multipattern)

    result = None
    for arg_val, arg_pattern in zip(tensors, patterns, strict=True):
        assert arg_val.dtype == torch.bool
        if result is None:
            result = rearrange(arg_val, arg_pattern)
        else:
            result = result & rearrange(arg_val, arg_pattern)
    return result


def und_self(mask: Tensor, pattern: str) -> Tensor:
    """
    Performs & (logical_and) for two replicas of the same tensor

    > und_self(mask, "b i, b j -> b 1 i j")
    is a better version of
    > und(mask, mask, "b i, b j -> b 1 i j")
    """
    return und(mask, mask, pattern)


# 255 is not an ASCII char
TENSORCODE_PAD_TOKEN = torch.iinfo(torch.uint8).max


@typecheck
def string_to_tensorcode(
    input: str,
    pad_to_length: int | None = None,
    device: torch.device | None = None,
) -> UInt8[Tensor, "l"]:
    """
    Converts an ASCII string to a tensor of integers.

    If pad_to_length is specified, the output tensor will have this length, and we add a
    special padding character if the tensor has less than the specified length.

    The minimum value of the output tensor is 0, and the maximum is 127 (excluding the
    padding token, which can be 255).
    """
    assert input.isascii(), "Expected input to be ASCII"
    ords = [ord(c) for c in input]

    tensorcode = torch.tensor(ords, dtype=torch.uint8, device=device)
    if pad_to_length is None:
        return tensorcode

    input_length = len(input)
    assert (
        pad_to_length >= input_length
    ), f"Expected {input_length=} to be shorter than {pad_to_length=} for {input=}"

    return F.pad(
        tensorcode,
        (0, pad_to_length - input_length),
        value=TENSORCODE_PAD_TOKEN,
    )


@typecheck
def tensorcode_to_string(tensor: UInt8[Tensor, "l"]) -> str:
    """
    Applies the inverse of the string_to_tensorcode function
    """
    assert tensor.device == torch.device("cpu")
    chars = [chr(i) for i in tensor if i != TENSORCODE_PAD_TOKEN]
    return "".join(chars)


@typecheck
def batch_tensorcode_to_string(
    tensor: UInt8[Tensor, "*dims l"],
) -> list[str]:
    tensor = rearrange(tensor, "... l -> (...) l")
    tensor = tensor[tensor.amax(dim=1) > 0, :]
    return [
        "".join(chr(i) for i in row if i != TENSORCODE_PAD_TOKEN)
        for row in tensor.tolist()
    ]


def unique_indexes(x: torch.Tensor, dim=-1, sorted: bool = True):
    """Implements return_index=True behavior for torch.unique.

    See https://numpy.org/doc/stable/reference/generated/numpy.unique.html for info and
    https://github.com/pytorch/pytorch/issues/36748 for context."""
    assert x.size(dim) > 0

    unique, inverse = torch.unique(x, return_inverse=True, sorted=True, dim=dim)
    perm = torch.arange(inverse.size(0), dtype=inverse.dtype, device=inverse.device)
    inverse, perm = inverse.flip([0]), perm.flip([0])
    inverse = inverse.new_empty(unique.size(dim)).scatter_(0, inverse, perm)
    if sorted:
        inverse = inverse.sort().values

    return unique, inverse


T = TypeVar("T")


# mypy is too angry when this function is directly annotated
def _move_data_to_device(x, device: torch.device):
    if x is None:
        return None
    if isinstance(x, (str, int, float, bool)):
        return x
    if isinstance(x, torch.Tensor):
        return x.to(device=device)
    elif isinstance(x, dict):
        return {k: move_data_to_device(v, device) for k, v in x.items()}
    elif isinstance(x, list):
        return [move_data_to_device(el, device) for el in x]
    elif isinstance(x, tuple):
        return tuple(move_data_to_device(el, device) for el in x)
    if is_dataclass(x):
        return replace(
            x,
            **{k: move_data_to_device(v, device) for k, v in asdict(x).items()},  # type: ignore
        )
    else:
        raise NotImplementedError(type(x))


def move_data_to_device(x: T, device: torch.device) -> T:
    return _move_data_to_device(x, device=device)


def set_seed(seed_sequence: list[int]) -> None:
    """
    Seeds numpy, torch, and Python.

    This function is heavily inspired by Lightning's pl_worker_init_function.
    """
    import random

    import numpy as np

    # Spawn distinct SeedSequences for the PyTorch PRNG and the stdlib random module
    np_ss = np.random.SeedSequence(seed_sequence)
    torch_ss, stdlib_ss = np_ss.spawn(2)

    # Seed numpy, use 128 bits (4 x 32-bit words)
    np.random.seed(np_ss.generate_state(4))

    # Seed torch
    torch.manual_seed(torch_ss.generate_state(1, dtype=np.uint64)[0])

    # Seed python, use 128 bits expressed as an integer
    stdlib_seed = (
        stdlib_ss.generate_state(2, dtype=np.uint64).astype(object) * [1 << 64, 1]
    ).sum()
    random.seed(stdlib_seed)

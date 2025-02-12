# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
"""Helper functions for aligning templates."""

import logging
from typing import Sequence

import torch
from torch import Tensor

from chai_lab.utils.typing import Shaped, typecheck

logger = logging.getLogger(__name__)


@typecheck
def align_1d(
    n_tokens: int,
    tensors: list[Tensor],
    indices: list[Tensor],
    pad_to_n_templates: int,
    default_dtype: torch.dtype,
    fill_value: int | float | bool = 0.0,
) -> Shaped[Tensor, "tmpl tok"]:
    """Align feature that is of shape (n_templates, n_tokens)."""
    n_temps = len(tensors)

    if pad_to_n_templates > 0:
        assert n_temps <= pad_to_n_templates
        n_temps = pad_to_n_templates

    retval = torch.full((n_temps, n_tokens), fill_value=fill_value, dtype=default_dtype)

    for i, (tensor, index) in enumerate(zip(tensors, indices, strict=True)):
        retval[i, index] = tensor

    return retval


@typecheck
def align_2d(
    n_tokens: int,
    tensors: list[Tensor],
    indices: list[Tensor],
    pad_to_n_templates: int,
    default_dtype: torch.dtype,
    addtl_dims: Sequence[int] = (),
) -> Shaped[Tensor, "tmpl tok tok *addtl"]:
    """Align feature that is of shape (n_templates, n_tokens, n_tokens)."""
    n_temps = len(tensors)
    assert n_temps == len(indices)

    if pad_to_n_templates > 0:
        assert n_temps <= pad_to_n_templates
        n_temps = pad_to_n_templates

    shape = (n_temps, n_tokens, n_tokens, *addtl_dims)
    retval = torch.zeros(shape, dtype=default_dtype)

    for i, (tensor, index) in enumerate(zip(tensors, indices, strict=True)):
        mask = torch.zeros(n_tokens, dtype=torch.bool)
        mask[index] = True
        a, b = torch.where(mask[:, None] * mask[None, :])
        retval[i, a, b] = tensor.flatten(0, 1)

    return retval

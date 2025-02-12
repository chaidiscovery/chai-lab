# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.utils.typing import Bool, typecheck


@typecheck
def subsample_msa_rows(
    mask: Bool[Tensor, "1 depth tokens"],
    select_n_rows: int = 4096,
    generator: torch.Generator | None = None,
) -> Bool[Tensor, "1 depth tokens"]:
    """Adjust masking to look at a random subset of msas.

    Returns input mask as-is if select_n_rows <= 0 or depth < select_n_rows."""
    nonnull_rows_mask = rearrange(mask.any(dim=-1), "1 d -> d")
    input_depth = nonnull_rows_mask.sum().item()
    if select_n_rows <= 0 or input_depth <= select_n_rows:
        return mask

    # Select from rows of the MSA that are not fully masked out
    (nonnull_row_indices,) = torch.where(nonnull_rows_mask)
    assert (n := nonnull_row_indices.numel()) > select_n_rows
    permuted = torch.randperm(n, device=mask.device, generator=generator)
    selected_row_indices = nonnull_row_indices[permuted[:select_n_rows]]

    # Create a mask for selected row indices
    selection_mask = torch.zeros_like(nonnull_rows_mask)
    selection_mask[selected_row_indices] = True
    selection_mask = repeat(selection_mask, "d -> 1 d 1")

    return mask & selection_mask

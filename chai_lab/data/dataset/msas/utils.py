# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange, reduce, repeat
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
    # Count the number of non-padding residues in each row of the MSA
    msa_sizes = rearrange(
        reduce(mask, "b depth tok -> b depth", reduction="sum"), "1 depth -> depth"
    )
    nonnull_rows_mask = msa_sizes > 0
    input_depth = nonnull_rows_mask.sum().item()
    if select_n_rows <= 0 or input_depth <= select_n_rows:
        return mask

    # Bias towards bigger hit MSAs; 0 size is automatically nulled out
    mask_ranking = msa_sizes * torch.rand(
        size=msa_sizes.shape,
        dtype=torch.float16,
        device=msa_sizes.device,
        generator=generator,
    )
    # Ascending sort -> choose the last (highest scoring) rows
    selected_row_indices = mask_ranking.argsort()[-select_n_rows:]

    # Create a mask for selected row indices
    selection_mask = torch.zeros_like(nonnull_rows_mask)
    selection_mask[selected_row_indices] = True

    return mask & repeat(selection_mask, "d -> 1 d 1")

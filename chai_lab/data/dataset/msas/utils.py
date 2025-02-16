# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging

import torch
from einops import rearrange, reduce
from torch import Tensor
from torch.nn import functional as F

from chai_lab.utils.typing import Bool, Float, typecheck


@typecheck
def _subsample_msa_rows(
    mask: Bool[Tensor, "1 depth tokens"],
    select_n_rows: int = 4096,
    generator: torch.Generator | None = None,
) -> Bool[Tensor, "depth"] | None:
    """Adjust masking to look at a random subset of msas.

    Returns None if select_n_rows <= 0 or depth < select_n_rows."""
    # Count the number of non-padding residues in each row of the MSA
    msa_sizes = rearrange(
        reduce(mask, "b depth tok -> b depth", reduction="sum"), "1 depth -> depth"
    )
    nonnull_rows_mask = msa_sizes > 0
    input_depth = nonnull_rows_mask.sum().item()
    if select_n_rows <= 0 or input_depth <= select_n_rows:
        return None

    # Bias towards bigger hit MSAs; 0 size is automatically nulled out
    mask_ranking = msa_sizes * torch.rand(
        size=msa_sizes.shape,
        dtype=torch.float16,
        device=msa_sizes.device,
        generator=generator,
    )
    # Ascending sort -> choose the last (highest scoring) rows
    selected_row_indices = mask_ranking.argsort()[-select_n_rows:]
    # We should never sample empty MSA rows
    assert not (~nonnull_rows_mask[selected_row_indices]).any()

    selection_mask = torch.zeros_like(nonnull_rows_mask)
    selection_mask[selected_row_indices] = True
    return selection_mask


@typecheck
def subsample_and_reorder_msa_feats_n_mask(
    feats: Float[Tensor, "1 depth tokens dim"],
    mask: Bool[Tensor, "1 depth tokens"],
    select_n_rows: int = 4096,
    generator: torch.Generator | None = None,
) -> tuple[Float[Tensor, "1 depth tokens dim"], Bool[Tensor, "1 depth tokens"]]:
    selection_mask = _subsample_msa_rows(
        mask=mask,
        select_n_rows=select_n_rows,
        generator=generator,
    )
    if selection_mask is None:  # No subsampling
        return feats, mask

    # Select the rows; where returns in order from top to bottom, preserving order
    (selection_idx,) = torch.where(selection_mask)
    logging.info(f"Subsampling {selection_idx.tolist()[:5]}...")
    (unselected_idx,) = torch.where(~selection_mask)
    combo_idx = torch.cat([selection_idx, unselected_idx])
    # Features are reordered, while mask is selected + padded
    feats_sampled = torch.index_select(feats, dim=1, index=combo_idx)
    mask_sampled = torch.index_select(mask, dim=1, index=selection_idx)
    # Every sampled row should have nonzero coverage
    assert mask_sampled.any(dim=-1).all()

    # Pad with zeros
    _, orig_depth, _ = mask.shape
    _, new_depth, _ = mask_sampled.shape
    assert (n_pad := orig_depth - new_depth) > 0
    # Padding is last dim, moving forward, e.g., for last two dimensions, it is:
    # (left, right, top, bottom)
    # [0, 0, 0, n_pad] ignores the token dim and pads out the depth dim
    return (
        feats_sampled,
        F.pad(mask_sampled, pad=[0, 0, 0, n_pad], value=False),
    )

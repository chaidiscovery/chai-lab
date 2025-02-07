# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

""""""

import logging
from pathlib import Path

import torch
from einops import reduce
from matplotlib import pyplot as plt
from torch import Tensor

from chai_lab.data import residue_constants as rc
from chai_lab.utils.typing import Int, UInt8, typecheck


@typecheck
def plot_msa(
    input_tokens: Int[Tensor, "n_tokens"],
    msa_tokens: UInt8[Tensor, "msa_depth n_tokens"],
    out_fname: Path,
    gap: str = "-",
    mask: str = ":",
    sort_by_identity: bool = True,
) -> Path:
    gap_idx = rc.residue_types_with_nucleotides.index(gap)
    mask_idx = rc.residue_types_with_nucleotides.index(mask)

    # Trim padding tokens (= pad in all alignments)
    token_is_pad = torch.all(msa_tokens == mask_idx, dim=0)
    msa_tokens = msa_tokens[:, ~token_is_pad]
    input_tokens = input_tokens[~token_is_pad]

    # Calculate sequence identity for each MSA sequence
    msa_seq_ident = (msa_tokens == input_tokens).float().mean(dim=-1)
    sort_idx = (
        torch.argsort(msa_seq_ident, descending=True)
        if sort_by_identity
        else torch.arange(msa_tokens.shape[0])
    )

    # Valid tokens are not padding and not a gap; we plot the valid tokens
    msa_tokens_is_valid = (msa_tokens != gap_idx) & (msa_tokens != mask_idx)
    msa_coverage = reduce(msa_tokens_is_valid.float(), "m t -> t", "mean")

    # Scale each of the MSA entries by its sequence identity for plotting
    msa_by_identity = msa_tokens_is_valid.float() * msa_seq_ident.unsqueeze(-1)
    msa_by_identity[~msa_tokens_is_valid] = torch.nan

    # Plotting
    fig, ax = plt.subplots(dpi=150)
    patch = ax.imshow(
        msa_by_identity[sort_idx],
        cmap="rainbow_r",
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )
    ax.set_aspect("auto")
    ax.set(ylabel="Sequences", xlabel="Positions")

    ax2 = ax.twinx()
    ax2.plot(msa_coverage, color="black")
    ax2.set(ylim=[0, 1], yticks=[])

    fig.colorbar(patch)
    fig.savefig(out_fname, bbox_inches="tight")
    logging.info(f"Saved MSA plot to {out_fname}")
    plt.close(fig)
    return out_fname

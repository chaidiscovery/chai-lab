""""""

import logging
from pathlib import Path
from typing import Any

import torch
from matplotlib import pyplot as plt
from torch import Tensor

from chai_lab.data import residue_constants as rc


def plot_msa(
    batch: dict[str, Any],
    out_fname: Path | str,
    gap: str = "-",
    mask: str = ":",
    sort_by_identity: bool = True,
):
    if "token_residue_type" not in batch or "msa_tokens" not in batch:
        logging.warning("No keys for MSA plotting; skipping...")
        return

    gap_idx = rc.residue_types_with_nucleotides.index(gap)
    mask_idx = rc.residue_types_with_nucleotides.index(mask)

    input_tokens = batch["token_residue_type"]  # Shape of (n_tokens,)
    assert isinstance(input_tokens, Tensor)
    assert input_tokens.ndim == 1
    n_tokens = input_tokens.numel()

    msa_tokens = batch["msa_tokens"]  # MSA depth, sequence length
    assert isinstance(msa_tokens, Tensor)
    assert (
        msa_tokens.ndim == 2 and msa_tokens.shape[-1] == n_tokens
    ), f"Bad shape: {msa_tokens.shape}"

    # Trim tokens that are all pad
    token_is_pad = torch.all(msa_tokens == mask_idx, dim=0)
    msa_tokens = msa_tokens[..., ~token_is_pad]
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
    msa_coverage = msa_tokens_is_valid.float().mean(dim=-2)

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

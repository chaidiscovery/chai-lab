# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import torch
from einops import rearrange, reduce, repeat
from torch import Tensor

from chai_lab.ranking.utils import expectation, get_chain_masks_and_asyms
from chai_lab.utils.tensor_utils import und
from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
@dataclass
class PTMScores:
    """
    complex_ptm: pTM score of the complex
    interface_ptm: ipTM score of the complex
    per_chain_ptm: pTM score for each chain in the complex
    per_chain_pair_iptm: ipTM score for each chain pair in the complex
    """

    complex_ptm: Float[Tensor, "..."]
    interface_ptm: Float[Tensor, "..."]
    per_chain_ptm: Float[Tensor, "... c"]
    per_chain_pair_iptm: Float[Tensor, "... c c"]


@typecheck
def tm_d0(n_tokens: Float[Tensor, "*dims"]) -> Float[Tensor, "*dims"]:
    """Compute TM-Score d0 from the number of tokens"""
    n_tokens = torch.clamp_min(n_tokens, 19)
    return 1.24 * (n_tokens - 15) ** (1.0 / 3) - 1.8


@typecheck
def _compute_ptm(
    logits: Float[Tensor, "... n n bins"],
    query_res_mask: Bool[Tensor, "... n"],
    query_has_frame_mask: Bool[Tensor, "... n"],
    key_res_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "bins"],
) -> Float[Tensor, "..."]:
    """
    Compute predicted TM score, normalized by the number of "key" tokens
    """
    num_key_tokens = reduce(key_res_mask, "... n -> ...", "sum").to(logits.dtype)
    # compute pairwise-TM normalized by the number of key tokens
    d0 = rearrange(tm_d0(num_key_tokens), "... -> ... 1")
    bin_weights: Float[Tensor, "bins"] = 1 / (1 + (bin_centers / d0) ** 2)
    # btm has shape (b,bins). Need to broadcast with probs
    # of shape (b,n,n,bins)
    bin_weights = rearrange(bin_weights, "... bins -> ... 1 1 bins")
    # determine key-query pairs with valid logits
    valid_pairs = und(
        query_has_frame_mask & query_res_mask, key_res_mask, "... i, ... j -> ... i j"
    )
    # compute per-pair expected TM scores
    expected_pair_tm = expectation(logits, bin_weights)
    # normalized scores by the number of key tokens
    num_key_tokens = rearrange(num_key_tokens, "... -> ... 1 1")
    qk_weights = valid_pairs.float() / torch.clamp_min(num_key_tokens, 1)
    # (b i j) -> (b i)
    query_key_tm = torch.sum(qk_weights * expected_pair_tm, dim=-1)
    # want to select the row with the most optimistic logits
    # and compute TM for this rows predicted alignment
    return torch.max(query_key_tm, dim=-1)[0]


@typecheck
def complex_ptm(
    pae_logits: Float[Tensor, "... n n n_bins"],
    token_exists_mask: Bool[Tensor, "... n"],
    valid_frames_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "n_bins"],
) -> Float[Tensor, "..."]:
    """Compute pTM score of the complex"""
    return _compute_ptm(
        logits=pae_logits,
        query_res_mask=token_exists_mask,
        query_has_frame_mask=valid_frames_mask,
        key_res_mask=token_exists_mask,
        bin_centers=bin_centers,
    )


@typecheck
def interface_ptm(
    pae_logits: Float[Tensor, "... n n n_bins"],
    token_exists_mask: Bool[Tensor, "... n"],
    valid_frames_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "n_bins"],
    token_asym_id: Int[Tensor, "... n"],
) -> Float[Tensor, "..."]:
    """Compute Interface pTM score

    ipTM is the max TM score over chains c \in C, restricting
    to interactions between c and C - {c}.
    """
    query_res_mask, _ = get_chain_masks_and_asyms(
        asym_id=token_asym_id, mask=token_exists_mask
    )

    per_chain_ptm = _compute_ptm(
        logits=rearrange(pae_logits, "... i j n_bins -> ... 1 i j n_bins"),
        query_res_mask=query_res_mask,
        query_has_frame_mask=rearrange(valid_frames_mask, "... n -> ... 1 n"),
        key_res_mask=~query_res_mask & rearrange(token_exists_mask, "... n -> ... 1 n"),
        bin_centers=bin_centers,
    )

    return torch.max(per_chain_ptm, dim=-1)[0]


@typecheck
def per_chain_pair_iptm(
    pae_logits: Float[Tensor, "... n n n_bins"],
    token_exists_mask: Bool[Tensor, "... n"],
    valid_frames_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "n_bins"],
    token_asym_id: Int[Tensor, "... n"],
    batched=False,
) -> tuple[Float[Tensor, "... n_chains n_chains"], Int[Tensor, "n_chains"]]:
    """Compute pairwise pTM score for each chain pair"""
    chain_mask, asyms = get_chain_masks_and_asyms(
        asym_id=token_asym_id, mask=token_exists_mask
    )
    c = asyms.numel()
    size = 32 * chain_mask.numel() ** 2 * c**2

    batched = batched and size < 2**32

    if not batched:
        # in the interest of saving memory we compute this in a for-loop
        results = []
        for i in range(c):
            result = _compute_ptm(
                logits=rearrange(pae_logits, "... i j n_bins -> ... 1 i j n_bins"),
                query_res_mask=repeat(chain_mask[..., i, :], "... n -> ... k n", k=c),
                query_has_frame_mask=rearrange(valid_frames_mask, "... n -> ... 1 n"),
                key_res_mask=chain_mask,
                bin_centers=bin_centers,
            )
            results.append(result)
        return torch.stack(results, dim=-2), asyms  # b, query_chain, key_chain
    else:
        # compute batched
        query_mask = repeat(chain_mask, "... c n -> ... c k n", k=c)
        key_mask = repeat(chain_mask, "... c n -> ... k c n", k=c)
        result = _compute_ptm(
            logits=rearrange(pae_logits, "... i j n_bins -> ... 1 1 i j n_bins"),
            query_res_mask=query_mask,
            query_has_frame_mask=rearrange(valid_frames_mask, "... n -> ... 1 1 n"),
            key_res_mask=key_mask,
            bin_centers=bin_centers,
        )
        return result, asyms


@typecheck
def per_chain_ptm(
    pae_logits: Float[Tensor, "... n n n_bins"],
    token_exists_mask: Bool[Tensor, "... n"],
    valid_frames_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "n_bins"],
    token_asym_id: Int[Tensor, "... n"],
) -> tuple[Float[Tensor, "... n_chains"], Int[Tensor, "n_chains"]]:
    """Computes pTM for each chain in the input"""
    chain_mask, unique_asyms = get_chain_masks_and_asyms(
        asym_id=token_asym_id, mask=token_exists_mask
    )
    per_chain_ptm = _compute_ptm(
        logits=rearrange(pae_logits, "... i j n_bins -> ... 1 i j n_bins"),
        query_res_mask=chain_mask,
        query_has_frame_mask=rearrange(valid_frames_mask, "... n -> ... 1 n"),
        key_res_mask=chain_mask,
        bin_centers=bin_centers,
    )
    return per_chain_ptm, unique_asyms


@typecheck
def get_scores(
    pae_logits: Float[Tensor, "... n n n_bins"],
    token_exists_mask: Bool[Tensor, "... n"],
    valid_frames_mask: Bool[Tensor, "... n"],
    bin_centers: Float[Tensor, "n_bins"],
    token_asym_id: Int[Tensor, "... n"],
) -> PTMScores:
    return PTMScores(
        complex_ptm=complex_ptm(
            pae_logits=pae_logits,
            token_exists_mask=token_exists_mask,
            valid_frames_mask=valid_frames_mask,
            bin_centers=bin_centers,
        ),
        interface_ptm=interface_ptm(
            pae_logits=pae_logits,
            token_exists_mask=token_exists_mask,
            valid_frames_mask=valid_frames_mask,
            bin_centers=bin_centers,
            token_asym_id=token_asym_id,
        ),
        per_chain_pair_iptm=per_chain_pair_iptm(
            pae_logits=pae_logits,
            token_exists_mask=token_exists_mask,
            valid_frames_mask=valid_frames_mask,
            bin_centers=bin_centers,
            token_asym_id=token_asym_id,
        )[0],
        per_chain_ptm=per_chain_ptm(
            pae_logits=pae_logits,
            token_exists_mask=token_exists_mask,
            valid_frames_mask=valid_frames_mask,
            bin_centers=bin_centers,
            token_asym_id=token_asym_id,
        )[0],
    )

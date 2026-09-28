# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import pytest
import torch

from chai_lab.ranking.rank import rank

_PAE_BINS = 32
_LDDT_BINS = 50


def _rank_inputs(chain_sizes: list[int]) -> dict:
    """Build minimal synthetic inputs for rank() with one atom per token.

    All tensors carry a leading batch dim of size 1.
    """
    n_tokens = sum(chain_sizes)
    token_asym_id = torch.cat(
        [
            torch.full((size,), asym_id, dtype=torch.int)
            for asym_id, size in enumerate(chain_sizes, start=1)
        ]
    )
    return {
        # Atoms spaced 5A apart so no clashes are possible
        "atom_coords": (
            torch.arange(n_tokens).unsqueeze(-1).repeat(1, 3).float() * 5.0
        ).unsqueeze(0),
        "atom_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "atom_token_index": torch.arange(n_tokens, dtype=torch.int).unsqueeze(0),
        "token_exists_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "token_asym_id": token_asym_id.unsqueeze(0),
        "token_entity_type": torch.zeros(1, n_tokens, dtype=torch.int),  # PROTEIN
        "token_valid_frames_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "lddt_logits": torch.zeros(1, n_tokens, _LDDT_BINS),
        "lddt_bin_centers": torch.linspace(0, 1, _LDDT_BINS),
        "pae_logits": torch.zeros(1, n_tokens, n_tokens, _PAE_BINS),
        "pae_bin_centers": torch.linspace(0, 31, _PAE_BINS),
    }


def test_single_chain_aggregate_score_is_ptm():
    """ipTM is undefined for a single chain; the aggregate must not blend it in (#327)."""
    ranking = rank(**_rank_inputs([10]))

    ptm = ranking.ptm_scores.complex_ptm
    assert ranking.ptm_scores.interface_ptm.item() == 0.0
    assert ranking.aggregate_score.item() == pytest.approx(ptm.item())
    # Guard the distinction: the buggy blend reports a much lower score
    assert ranking.aggregate_score.item() != pytest.approx(0.2 * ptm.item())


def test_multi_chain_aggregate_score_unchanged():
    """Multi-chain inputs keep the 0.2*pTM + 0.8*ipTM blend."""
    ranking = rank(**_rank_inputs([6, 6]))

    expected = (
        0.2 * ranking.ptm_scores.complex_ptm.item()
        + 0.8 * ranking.ptm_scores.interface_ptm.item()
    )
    assert ranking.ptm_scores.interface_ptm.item() > 0.0
    assert ranking.aggregate_score.item() == pytest.approx(expected)

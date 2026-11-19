# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import numpy as np
import torch

from chai_lab.chai1 import write_scores_npz
from chai_lab.ranking.rank import rank

_PAE_BINS = 32
_LDDT_BINS = 50


def _rank_inputs(n_tokens: int) -> dict:
    """Minimal synthetic single-chain inputs for rank(), one atom per token."""
    return {
        # Atoms spaced 5A apart so no clashes are possible
        "atom_coords": (
            torch.arange(n_tokens).unsqueeze(-1).repeat(1, 3).float() * 5.0
        ).unsqueeze(0),
        "atom_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "atom_token_index": torch.arange(n_tokens, dtype=torch.int).unsqueeze(0),
        "token_exists_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "token_asym_id": torch.ones(1, n_tokens, dtype=torch.int),
        "token_entity_type": torch.zeros(1, n_tokens, dtype=torch.int),  # PROTEIN
        "token_valid_frames_mask": torch.ones(1, n_tokens, dtype=torch.bool),
        "lddt_logits": torch.zeros(1, n_tokens, _LDDT_BINS),
        "lddt_bin_centers": torch.linspace(0, 1, _LDDT_BINS),
        "pae_logits": torch.zeros(1, n_tokens, n_tokens, _PAE_BINS),
        "pae_bin_centers": torch.linspace(0, 31, _PAE_BINS),
    }


def test_scores_npz_contains_confidence_matrices(tmp_path):
    """#348: Python-mode runs must persist PAE/PDE/pLDDT, not just scalars."""
    n_tokens = 8
    n_atoms = 8
    ranking = rank(**_rank_inputs(n_tokens))

    out = tmp_path / "scores.model_idx_0.npz"
    write_scores_npz(
        out,
        ranking,
        pae_scores=torch.rand(1, n_tokens, n_tokens),
        pde_scores=torch.rand(1, n_tokens, n_tokens),
        plddt_scores_atom=torch.rand(1, n_atoms),
    )

    loaded = np.load(out)
    assert loaded["pae"].shape == (1, n_tokens, n_tokens)
    assert loaded["pde"].shape == (1, n_tokens, n_tokens)
    assert loaded["plddt"].shape == (1, n_atoms)
    # scalar metrics still present
    assert "aggregate_score" in loaded
    assert "ptm" in loaded

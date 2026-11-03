# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Callable

import pytest
import torch

from chai_lab.ranking.ptm import get_scores, tm_d0

SCORE_FIELDS = [
    "complex_ptm",
    "interface_ptm",
    "per_chain_ptm",
    "per_chain_pair_iptm",
]


def _mismatch_msg(field: str, detail: str) -> Callable[[str], str]:
    def build(generated: str) -> str:
        return f"{field} changed for {detail}:\n{generated}"

    return build


def _uniform_logit_inputs(chain_sizes: list[int], n_bins: int = 16) -> dict:
    """Inputs whose softmax is exactly representable in every float dtype.

    Constant logits give a uniform distribution of 1/n_bins; with n_bins a power
    of two that value is exact in bf16, fp16 and float32 alike. So the scores
    computed from narrow logits and from float32 logits must agree exactly, and
    any difference is attributable to how the token counts are handled rather
    than to the precision of the probabilities themselves.
    """
    n_tokens = sum(chain_sizes)
    return dict(
        pae_logits=torch.zeros(1, n_tokens, n_tokens, n_bins),
        token_asym_id=torch.cat(
            [torch.full((s,), i + 1) for i, s in enumerate(chain_sizes)]
        )
        .int()
        .unsqueeze(0),
        token_exists_mask=torch.ones(1, n_tokens, dtype=torch.bool),
        valid_frames_mask=torch.ones(1, n_tokens, dtype=torch.bool),
        bin_centers=torch.linspace(0.0, 32.0, n_bins),
    )


def _random_inputs(seed: int, chain_sizes: list[int], n_bins: int = 16) -> dict:
    g = torch.Generator().manual_seed(seed)
    n_tokens = sum(chain_sizes)
    return dict(
        pae_logits=torch.randn(1, n_tokens, n_tokens, n_bins, generator=g),
        token_asym_id=torch.cat(
            [torch.full((s,), i + 1) for i, s in enumerate(chain_sizes)]
        )
        .int()
        .unsqueeze(0),
        token_exists_mask=torch.ones(1, n_tokens, dtype=torch.bool),
        valid_frames_mask=torch.rand(1, n_tokens, generator=g) > 0.2,
        bin_centers=torch.linspace(0.0, 32.0, n_bins),
    )


def test_tm_d0_is_not_evaluated_in_low_precision():
    """d0 must be computed in at least float32.

    tm_d0 is 1.24 * (n - 15)**(1/3) - 1.8. Evaluated in bfloat16 -- which has an
    8-bit mantissa -- this carries a relative error of several 1e-3 even for
    token counts that are themselves exactly representable.
    """
    for n in [128, 256, 384, 512, 1024]:
        narrow = tm_d0(torch.tensor(float(n), dtype=torch.bfloat16)).float()
        wide = tm_d0(torch.tensor(float(n), dtype=torch.float32))
        # Confirms the premise: bf16 really is too coarse to evaluate d0 in
        assert not torch.isclose(narrow, wide, atol=1e-4, rtol=0), (
            f"bf16 d0 for n={n} is unexpectedly accurate; this test no longer "
            "demonstrates why counts must be widened"
        )


@pytest.mark.parametrize("dtype", [torch.bfloat16, torch.float16])
@pytest.mark.parametrize(
    "chain_sizes",
    [
        [256, 128],  # both sizes exactly representable in bf16
        [275, 40],  # 275 -> 276 in bf16
        [257, 40],  # 257 -> 256 in bf16
        [431, 40],  # 431 -> 432 in bf16
    ],
)
def test_scores_do_not_depend_on_logits_dtype(
    dtype: torch.dtype, chain_sizes: list[int]
):
    """Scores must not change when identical probabilities arrive as bf16/fp16.

    The token counts feed d0 and the 1/|K| normalization. Deriving them from the
    logits dtype makes both the counts and the d0 arithmetic low-precision, so
    the scores shift by ~1e-3 for narrow logits even though the underlying
    distribution is bit-identical.
    """
    kwargs = _uniform_logit_inputs(chain_sizes)
    narrow = dict(kwargs, pae_logits=kwargs["pae_logits"].to(dtype))

    # Control: the probabilities really are identical in the two dtypes, so any
    # difference in the scores cannot come from the softmax
    torch.testing.assert_close(
        torch.softmax(narrow["pae_logits"], dim=-1).float(),
        torch.softmax(kwargs["pae_logits"], dim=-1),
        atol=0,
        rtol=0,
    )

    from_narrow = get_scores(**narrow)
    from_float32 = get_scores(**kwargs)
    for field in SCORE_FIELDS:
        torch.testing.assert_close(
            getattr(from_narrow, field).float(),
            getattr(from_float32, field).float(),
            atol=1e-6,
            rtol=0,
            msg=_mismatch_msg(field, f"{dtype} logits with chain sizes {chain_sizes}"),
        )


@pytest.mark.parametrize("chain_sizes", [[20], [10, 14], [8, 10, 12], [6] * 5])
def test_score_shapes_and_ranges(chain_sizes: list[int]):
    scores = get_scores(
        **_random_inputs(seed=len(chain_sizes), chain_sizes=chain_sizes)
    )
    n_chains = len(chain_sizes)

    assert scores.per_chain_ptm.shape == (1, n_chains)
    assert scores.per_chain_pair_iptm.shape == (1, n_chains, n_chains)

    for field in SCORE_FIELDS:
        value = getattr(scores, field)
        assert torch.all(value >= 0) and torch.all(value <= 1), (
            f"{field} outside [0, 1]: {value}"
        )


def test_single_chain_has_no_interface():
    """With one chain there is nothing to interface with, so ipTM is zero."""
    scores = get_scores(**_random_inputs(seed=3, chain_sizes=[20]))
    torch.testing.assert_close(scores.interface_ptm, torch.zeros(1))


def test_per_chain_ptm_matches_self_pair():
    """per-chain pTM is the query==key entry of the chain-pair matrix."""
    scores = get_scores(**_random_inputs(seed=4, chain_sizes=[8, 10, 12]))
    torch.testing.assert_close(
        scores.per_chain_ptm,
        torch.diagonal(scores.per_chain_pair_iptm, dim1=-2, dim2=-1),
    )

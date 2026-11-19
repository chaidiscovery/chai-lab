# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


def get_sep_bins(max_offset: int) -> list[float]:
    bins = torch.arange(-max_offset, max_offset + 1).float()
    return bins.tolist()


SMALL_SEP_BINS = get_sep_bins(32)


class RelativeSequenceSeparation(FeatureGenerator):
    def __init__(
        self,
        sep_bins: list[int] | list[float] | None = None,
        num_bins: int | None = None,
    ):
        """Relative Sequence Separation Encoding"""
        sep_bins = get_sep_bins(num_bins) if num_bins is not None else sep_bins
        sep_bins = sep_bins if sep_bins is not None else SMALL_SEP_BINS
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            num_classes=len(sep_bins) + 2,
            can_mask=False,
        )
        self.sep_bins = torch.tensor(sep_bins)

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            residue_index=batch["inputs"]["token_residue_index"].long(),
            asym_id=batch["inputs"]["token_asym_id"].long(),
            cyclic_period=batch["inputs"]["token_cyclic_period"].long(),
        )

    @typecheck
    def _generate(
        self,
        residue_index: Int[Tensor, "b n"],
        asym_id: Int[Tensor, "b n"],
        cyclic_period: Int[Tensor, "b n"],
    ) -> Tensor:
        rel_sep, rel_chain = map(
            lambda x: rearrange(x, "b n -> b n 1") - rearrange(x, "b n -> b 1 n"),
            (residue_index, asym_id),
        )
        same_chain_mask = rel_chain == 0
        if torch.any(cyclic_period > 0):
            rel_period = torch.maximum(
                rearrange(cyclic_period, "b n -> b n 1"),
                rearrange(cyclic_period, "b n -> b 1 n"),
            )
            safe_period = torch.where(rel_period > 0, rel_period, torch.ones_like(rel_period))
            cyclic_pair_mask = same_chain_mask & (rel_period > 0)
            wrapped_rel_sep = (
                rel_sep - rel_period * torch.round(rel_sep.float() / safe_period.float())
            ).long()
            rel_sep = torch.where(cyclic_pair_mask, wrapped_rel_sep, rel_sep)

        encoded_feat = torch.searchsorted(
            self.sep_bins.to(rel_sep.device),
            rel_sep + 1e-4,  # add small epsilon bc. bins are chosen by leftmost index
        )
        # mask inter-chain sep
        encoded_feat[~same_chain_mask] = self.num_classes - 1
        return self.make_feature(encoded_feat.unsqueeze(-1))

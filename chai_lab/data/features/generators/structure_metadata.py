# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import repeat
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.defaults import default
from chai_lab.utils.typing import Bool, Float, typecheck

DEFAULT_BFACTOR_BINS = [140.0]

DEFAULT_PLDDT_BINS = [0.3, 0.7]


class IsDistillation(FeatureGenerator):
    def __init__(self):
        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            is_distillation=batch["inputs"]["is_distillation"],
            token_exists_mask=batch["inputs"]["token_exists_mask"],
        )

    @typecheck
    def _generate(
        self,
        is_distillation: Bool[Tensor, "b 1"],
        token_exists_mask: Bool[Tensor, "b n"],
    ) -> Tensor:
        _, n = token_exists_mask.shape
        is_distillation = repeat(is_distillation, "b 1 -> b n 1", n=n).to(torch.uint8)
        return self.make_feature(data=is_distillation)


class TokenBFactor(FeatureGenerator):
    def __init__(
        self,
        include_prob: float = 1.0,
        bins: list[float] | None = None,
    ):
        self.bins = torch.tensor(default(bins, DEFAULT_BFACTOR_BINS))

        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=len(self.bins) + 1,
            mult=1,
        )
        self.include_prob = include_prob

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            token_b_factor=batch["inputs"]["token_b_factor_or_plddt"],
            is_distillation=batch["inputs"]["is_distillation"],
            token_exists_mask=batch["inputs"]["token_exists_mask"],
        )

    @typecheck
    def _generate(
        self,
        token_b_factor: Float[Tensor, "b n"],
        is_distillation: Bool[Tensor, "b 1"],
        token_exists_mask: Bool[Tensor, "b n"],
    ) -> Tensor:
        _, n = token_exists_mask.shape

        include_mask = (
            torch.rand_like(is_distillation, dtype=torch.float) <= self.include_prob
        )

        # this feature is not defined for distillation data
        mask = (
            repeat(~is_distillation, "b 1 -> b n", n=n)
            & token_exists_mask
            & repeat(include_mask, "b 1 -> b n", n=n)
        )

        feat = torch.searchsorted(self.bins.to(is_distillation.device), token_b_factor)
        feat.masked_fill_(~mask, self.mask_value)

        return self.make_feature(data=feat.unsqueeze(-1))


class TokenPLDDT(FeatureGenerator):
    def __init__(
        self,
        include_prob: float = 1.0,
        bins: list[float] | None = None,
    ):
        self.bins = torch.tensor(default(bins, DEFAULT_PLDDT_BINS))

        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=len(self.bins) + 1,
            mult=1,
        )
        self.include_prob = include_prob

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            token_plddt=batch["inputs"]["token_b_factor_or_plddt"],
            is_distillation=batch["inputs"]["is_distillation"],
            token_exists_mask=batch["inputs"]["token_exists_mask"],
        )

    @typecheck
    def _generate(
        self,
        token_plddt: Float[Tensor, "b n"],
        is_distillation: Bool[Tensor, "b 1"],
        token_exists_mask: Bool[Tensor, "b n"],
    ) -> Tensor:
        _, n = token_exists_mask.shape

        include_mask = (
            torch.rand_like(is_distillation, dtype=torch.float) <= self.include_prob
        )

        # this feature is defined ONLY for distillation data
        mask = (
            repeat(is_distillation, "b 1 -> b n", n=n)
            & token_exists_mask
            & repeat(include_mask, "b 1 -> b n", n=n)
        )

        feat = torch.searchsorted(self.bins.to(is_distillation.device), token_plddt)
        feat.masked_fill_(~mask, self.mask_value)

        return self.make_feature(data=feat.unsqueeze(-1))

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class RelativeTokenSeparation(FeatureGenerator):
    def __init__(
        self,
        # using 16 for default here since values beyond this are very rare.
        r_max: int = 16,
    ):
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            num_classes=2 * r_max + 3,
            can_mask=False,
        )
        self.r_max = r_max

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            token_index=batch["inputs"]["token_index"],
            token_residue_index=batch["inputs"]["token_residue_index"],
            token_asym_id=batch["inputs"]["token_asym_id"],
        )

    @typecheck
    def _generate(
        self,
        token_index: Int[Tensor, "b n"],
        token_residue_index: Int[Tensor, "b n"],
        token_asym_id: Int[Tensor, "b n"],
    ) -> Tensor:
        rel_sep, rel_residue, rel_chain = map(
            lambda x: rearrange(x, "b n -> b n 1") - rearrange(x, "b n -> b 1 n"),
            (token_index, token_residue_index, token_asym_id),
        )

        mask = (rel_residue == 0) & (rel_chain == 0)

        rel_sep = torch.clamp(rel_sep + self.r_max, 0, 2 * self.r_max + 1)
        # zero inter-residue and inter-chain
        rel_sep = rel_sep.masked_fill(~mask, 2 * self.r_max + 2)

        return self.make_feature(rel_sep.unsqueeze(-1))

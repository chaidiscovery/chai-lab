# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class ChainIsCropped(FeatureGenerator):
    def __init__(
        self,
    ):
        """Chain-level feature that indicates if a chain has been cropped"""
        super().__init__(
            ty=FeatureType.TOKEN,
            can_mask=False,
            encoding_ty=EncodingType.IDENTITY,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            token_asym_id=batch["inputs"]["token_asym_id"].long(),
        )

    @typecheck
    def _generate(
        self,
        token_asym_id: Int[Tensor, "b n"],
    ) -> Tensor:
        return self.make_feature(torch.zeros_like(token_asym_id).unsqueeze(-1))

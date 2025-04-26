# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class ResidueType(FeatureGenerator):
    def __init__(
        self,
        num_res_ty: int = 22,  # 20AA + gap + X
        key: str = "aatype",
    ):
        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=num_res_ty,
            mult=1,
        )
        self.key = key

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(aatype=batch["inputs"][self.key].long())

    @typecheck
    def _generate(self, aatype: Int[Tensor, "b n"]) -> Tensor:
        """see super class"""
        seq_emb = aatype.clone()
        return self.make_feature(data=seq_emb.unsqueeze(-1))

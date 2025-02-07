# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class AtomElementOneHot(FeatureGenerator):
    def __init__(
        self,
        max_atomic_num: int = 128,
    ):
        super().__init__(
            ty=FeatureType.ATOM,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=max_atomic_num + 1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(atomic_numbers=batch["inputs"]["atom_ref_element"])

    @typecheck
    def _generate(self, atomic_numbers: Int[Tensor, "b n"]) -> Tensor:
        """see super class"""
        return self.make_feature(
            data=torch.clamp(atomic_numbers, max=self.num_classes).unsqueeze(-1),
        )

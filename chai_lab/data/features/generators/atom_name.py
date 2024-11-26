# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class AtomNameOneHot(FeatureGenerator):
    def __init__(
        self,
        num_chars: int = 64,
    ):
        super().__init__(
            ty=FeatureType.ATOM,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=num_chars,
            mult=4,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(atom_name_chars=batch["inputs"]["atom_ref_name_chars"])

    @typecheck
    def _generate(self, atom_name_chars: Int[Tensor, "b n 4"]) -> Tensor:
        """see super class"""
        return self.make_feature(data=atom_name_chars)

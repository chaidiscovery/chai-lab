# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import numpy as np
import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Bool, Int, typecheck


class ResidueType(FeatureGenerator):
    def __init__(
        self,
        min_corrupt_prob: float = 0.0,
        max_corrupt_prob: float = 0.0,
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
        self.min_corrupt_prob = min_corrupt_prob
        self.max_corrupt_prob = max_corrupt_prob
        self.key = key

    @typecheck
    def _corrupt_seq(
        self, sequence: Int[Tensor, "... n"]
    ) -> tuple[Int[Tensor, "... n"], Bool[Tensor, "... n"]]:
        """Corrupt the sequence with the given probability"""
        corrupt_prob = np.random.uniform(
            low=self.min_corrupt_prob, high=self.max_corrupt_prob
        )
        corrupt_mask = torch.rand_like(sequence.float()) < corrupt_prob
        corrupt_aas = torch.randint_like(
            corrupt_mask[corrupt_mask].long(), high=self.num_classes - 1
        )
        corrupt_sequence = sequence.clone()
        corrupt_sequence[corrupt_mask] = corrupt_aas
        return corrupt_sequence, corrupt_mask

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(aatype=batch["inputs"][self.key].long())

    @typecheck
    def _generate(self, aatype: Int[Tensor, "b n"]) -> Tensor:
        """see super class"""
        seq_emb = aatype.clone()
        seq_emb, _corrupt_mask = self._corrupt_seq(seq_emb)
        return self.make_feature(data=seq_emb.unsqueeze(-1))

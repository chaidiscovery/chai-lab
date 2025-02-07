# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator

N_COORDS = 3  # 3 coords: x, y, z


class RefPos(FeatureGenerator):
    """Provides reference position of atom"""

    def __init__(self):
        super().__init__(
            ty=FeatureType.ATOM,
            encoding_ty=EncodingType.IDENTITY,
            mult=1,
            num_classes=N_COORDS,
            can_mask=False,  # we expect to always have valid pos?
        )

    def generate(self, batch: dict) -> torch.Tensor:
        original_pos = batch["inputs"]["atom_ref_pos"]
        feat = original_pos / 10.0  # better scale for embedding
        assert torch.amax(feat.norm(dim=-1)) < 100.0, "wrong scale!"
        assert feat.ndim == 3
        assert feat.shape[-1] == N_COORDS

        return self.make_feature(data=feat)

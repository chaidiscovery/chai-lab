# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

import chai_lab.data.features.feature_utils as futils
from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator


class Identity(FeatureGenerator):
    def __init__(
        self,
        key: str,
        ty: FeatureType,
        dim: int,
        can_mask: bool = True,
    ):
        super().__init__(
            ty=ty,
            encoding_ty=EncodingType.IDENTITY,
            mult=1,
            num_classes=dim,
            can_mask=can_mask,
        )
        self.key = key
        self.dim = dim

    def generate(self, batch: dict) -> Tensor:
        feat = futils.get_entry_for_key(batch, self.key)

        if feat.ndim == 2:  # scalar feature
            assert self.dim == 1
            feat = rearrange(feat, "b n -> b n 1")
        elif feat.ndim == 3:
            # feature made from sequence-wise vectors (shape b,n,d)
            assert self.dim == feat.shape[-1]
        else:
            raise ValueError(
                f"Input to feature generator has ndim={feat.ndim}, shape {feat.shape}"
            )

        if self.can_mask:  # append position for mask token if feat can be masked
            feat = torch.cat((feat, torch.zeros_like(feat)[..., :1]), dim=-1)
        return self.make_feature(data=feat)

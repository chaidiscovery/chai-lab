# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class RelativeEntity(FeatureGenerator):
    def __init__(self):
        """Relative Entity Encoding

        See algorithm 5 of AF-Multimer
        """
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            num_classes=3,
            can_mask=False,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            entity_id=batch["inputs"]["token_entity_id"].long(),
        )

    @typecheck
    def _generate(
        self,
        entity_id: Int[Tensor, "b n"],
    ) -> Tensor:
        # remap unique sym_id values to 0,n-1
        _, entity_id_from_zero = torch.unique(
            entity_id, sorted=True, return_inverse=True
        )

        rel_entity = rearrange(entity_id_from_zero, "b n -> b n 1") - rearrange(
            entity_id_from_zero, "b n -> b 1 n"
        )
        rel_entity = torch.clamp(rel_entity + 1, 0, 2)
        assert rel_entity.dtype == torch.long
        return self.make_feature(rel_entity.unsqueeze(-1))

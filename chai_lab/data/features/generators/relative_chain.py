# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Int, typecheck


class RelativeChain(FeatureGenerator):
    def __init__(
        self,
        s_max: int = 2,
    ):
        """Relative Entity Encoding

        See algorithm 5 of AF-Multimer
        """
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            num_classes=2 * s_max + 2,
            can_mask=False,
        )
        self.s_max = s_max

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            entity_id=batch["inputs"]["token_entity_id"].long(),
            sym_id=batch["inputs"]["token_sym_id"].long(),
        )

    @typecheck
    def _generate(
        self,
        entity_id: Int[Tensor, "b n"],
        sym_id: Int[Tensor, "b n"],
    ) -> Tensor:
        # remap unique sym_id values to 0,n-1
        _, sym_ids_from_zero = torch.unique(sym_id, sorted=True, return_inverse=True)

        rel_entity, rel_chain = map(
            lambda x: rearrange(x, "b n -> b n 1") - rearrange(x, "b n -> b 1 n"),
            (entity_id, sym_ids_from_zero),
        )
        # within an entity, determine relative chain
        rel_chain = torch.clamp(rel_chain + self.s_max, 0, 2 * self.s_max)
        # mask out inter-entity features
        rel_chain[rel_entity != 0] = 2 * self.s_max + 1
        return self.make_feature(rel_chain.unsqueeze(-1))

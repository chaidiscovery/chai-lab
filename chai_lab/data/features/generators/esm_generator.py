# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.typing import Float, typecheck


class ESMEmbeddings(FeatureGenerator):
    def __init__(
        self,
        ty: FeatureType = FeatureType.TOKEN,
    ):
        super().__init__(
            ty=ty,
            encoding_ty=EncodingType.ESM,
            can_mask=False,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            esm_embeddings=batch["inputs"]["esm_embeddings"],
        )

    @typecheck
    def _generate(
        self,
        esm_embeddings: Float[Tensor, "batch num_tokens d_emb"],
    ) -> Tensor:
        return self.make_feature(data=esm_embeddings)

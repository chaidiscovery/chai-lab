# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Any

import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.features.token_utils import get_centre_positions_and_mask
from chai_lab.utils.tensor_utils import cdist
from chai_lab.utils.typing import Bool, Float, Int, typecheck


class TokenCenterDistance(FeatureGenerator):
    def __init__(
        self,
        dist_bins: list[float] | None = None,
    ):
        dist_bins = dist_bins if dist_bins is not None else [0.0, 4.0, 8.0, 12.0, 16.0]
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            # one of dist_bins of rbf_radii is not None.
            num_classes=len(dist_bins) + 1,
            mult=1,
            can_mask=True,
        )

        # maintain consistent orders
        self.dist_bins = torch.tensor(dist_bins)

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            all_atom_positions=batch["inputs"]["atom_gt_coords"],
            all_atom_mask=batch["inputs"]["atom_exists_mask"],
            token_single_mask=batch["inputs"]["token_exists_mask"],
            token_center_atom_index=batch["inputs"]["token_centre_atom_index"].long(),
        )

    @typecheck
    def _generate(
        self,
        all_atom_positions=Float[Tensor, "b a 3"],
        all_atom_mask=Bool[Tensor, "b a"],
        token_single_mask=Bool[Tensor, "b n"],
        token_center_atom_index=Int[Tensor, "b n"],
    ) -> Tensor:
        """see super class"""
        center_atom_coords, center_atom_mask = get_centre_positions_and_mask(
            atom_gt_coords=all_atom_positions,
            atom_exists_mask=all_atom_mask,
            token_centre_atom_index=token_center_atom_index,
            token_exists_mask=token_single_mask,
        )
        feat = torch.searchsorted(
            self.dist_bins.to(center_atom_coords.device), cdist(center_atom_coords)
        )
        center_atom_pair_exists = torch.einsum(
            "b i, b j -> b i j", center_atom_mask, center_atom_mask
        )
        feat.masked_fill_(~center_atom_pair_exists, self.mask_value)
        return self.make_feature(feat.unsqueeze(-1))

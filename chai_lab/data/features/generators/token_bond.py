# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Any

import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.tensor_utils import und_self
from chai_lab.utils.typing import Bool, Int, typecheck


class TokenBondRestraint(FeatureGenerator):
    def __init__(self):
        """Generates features for covalent bonds between atoms."""
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            token_exists_mask=batch["inputs"]["token_exists_mask"],
            atom_token_index=batch["inputs"]["atom_token_index"].long(),
            atom_covalent_bond_indices=batch["inputs"]["atom_covalent_bond_indices"],
        )

    @typecheck
    def apply_mask(self, feature: Tensor, mask: Tensor, mask_ty: FeatureType) -> Tensor:
        # override masking behavior - just return the unmasked feature
        return feature

    @typecheck
    def _generate(
        self,
        token_exists_mask: Bool[Tensor, "b n"],
        atom_token_index: Int[Tensor, "b a"],
        atom_covalent_bond_indices: list[
            tuple[Int[Tensor, "bonds"], Int[Tensor, "bonds"]]
        ],
    ) -> Tensor:
        token_pair_mask = und_self(token_exists_mask, "b i, b j -> b i j")
        bond_feature = torch.zeros_like(token_pair_mask.float())

        for batch_idx, (left_indices, right_indices) in enumerate(
            atom_covalent_bond_indices
        ):
            # convert from atom index to token index
            left_token_indices = torch.gather(
                atom_token_index[batch_idx], dim=0, index=left_indices
            )
            right_token_indices = torch.gather(
                atom_token_index[batch_idx], dim=0, index=right_indices
            )
            bond_feature[batch_idx][left_token_indices, right_token_indices] = 1

        return self.make_feature(bond_feature.unsqueeze(-1))

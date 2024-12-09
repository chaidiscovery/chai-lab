# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Any, Literal

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.tensor_utils import cdist
from chai_lab.utils.typing import Bool, Float, Int, typecheck

_VALID_ENCODING_TYPES = [
    EncodingType.IDENTITY,
]
DEFAULT_ONE_HOT_DIST_BINS = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0]
DEFAULT_RBF_DIST_BINS = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]


class BlockedAtomPairDistances(FeatureGenerator):
    transform: Literal["none", "inverse_squared"]

    def __init__(
        self,
        encoding_ty: EncodingType = EncodingType.IDENTITY,
        transform: Literal["none", "inverse_squared"] = "inverse_squared",
    ):
        assert (
            encoding_ty in _VALID_ENCODING_TYPES
        ), f"invalid encoding type: {encoding_ty}"

        # initialize superclass after augmenting input params =O.
        super().__init__(
            ty=FeatureType.ATOM_PAIR,
            encoding_ty=encoding_ty,
            # one of dist_bins of rbf_radii is not None.
            num_classes=1,
            mult=1,
            can_mask=True,
        )
        self.transform = transform

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            atom_ref_pos=batch["inputs"]["atom_ref_pos"],
            atom_ref_mask=batch["inputs"]["atom_ref_mask"],
            atom_ref_space_uid=batch["inputs"]["atom_ref_space_uid"],
            q_idces=batch["inputs"]["block_atom_pair_q_idces"],
            kv_idces=batch["inputs"]["block_atom_pair_kv_idces"],
            block_atom_pair_mask=batch["inputs"]["block_atom_pair_mask"],
        )

    @typecheck
    def _generate(
        self,
        atom_ref_pos: Float[Tensor, "b n 3"],
        atom_ref_mask: Bool[Tensor, "b n"],
        atom_ref_space_uid: Int[Tensor, "b n"],
        q_idces: Int[Tensor, "bl bl_q"],
        kv_idces: Int[Tensor, "bl bl_kv"],
        block_atom_pair_mask: Bool[Tensor, "b bl bl_q bl_kv"],
    ) -> Tensor:
        """see super class"""

        blocked_feat, blocked_mask = get_blocked_atom_pair_dists(
            atom_ref_pos,
            atom_ref_space_uid,
            q_idces,
            kv_idces,
            block_atom_pair_mask,
        )

        if self.transform == "inverse_squared":
            blocked_feat = 1 / (1 + blocked_feat**2)

        # return (B, n, n, 2) where ...,0 is the feature
        # and ...,1 indicates if the value is masked
        # because 0.0 has a meaning as a distance

        blocked_feat = blocked_feat.unsqueeze(-1)
        blocked_mask = blocked_mask.unsqueeze(-1).float()

        return self.make_feature(
            torch.cat(
                [blocked_feat, blocked_mask],
                dim=-1,
            )
        )


class BlockedAtomPairDistogram(FeatureGenerator):
    dist_bins: Tensor

    def __init__(
        self,
        dist_bins: list[float] | None = None,
        encoding_ty: EncodingType = EncodingType.ONE_HOT,
    ):
        if dist_bins is None and encoding_ty == EncodingType.ONE_HOT:
            dist_bins = DEFAULT_ONE_HOT_DIST_BINS
        elif dist_bins is None and encoding_ty == EncodingType.RBF:
            dist_bins = DEFAULT_RBF_DIST_BINS
        assert dist_bins is not None, "must provide dist_bins"

        # initialize superclass after augmenting input params =O.
        super().__init__(
            ty=FeatureType.ATOM_PAIR,
            encoding_ty=encoding_ty,
            # one of dist_bins of rbf_radii is not None.
            num_classes=len(dist_bins) + 1,
            mult=1,
            can_mask=True,
        )
        self.dist_bins = torch.tensor(dist_bins)

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            atom_ref_pos=batch["inputs"]["atom_ref_pos"],
            atom_ref_mask=batch["inputs"]["atom_ref_mask"],
            atom_ref_space_uid=batch["inputs"]["atom_ref_space_uid"],
            q_idces=batch["inputs"]["block_atom_pair_q_idces"],
            kv_idces=batch["inputs"]["block_atom_pair_kv_idces"],
            block_atom_pair_mask=batch["inputs"]["block_atom_pair_mask"],
        )

    @typecheck
    def _generate(
        self,
        atom_ref_pos: Float[Tensor, "b n 3"],
        atom_ref_mask: Bool[Tensor, "b n"],
        atom_ref_space_uid: Int[Tensor, "b n"],
        q_idces: Int[Tensor, "bl bl_q"],
        kv_idces: Int[Tensor, "bl bl_kv"],
        block_atom_pair_mask: Bool[Tensor, "b bl bl_q bl_kv"],
    ) -> Tensor:
        """see super class"""
        feat, mask = get_blocked_atom_pair_dists(
            atom_ref_pos,
            atom_ref_space_uid,
            q_idces,
            kv_idces,
            block_atom_pair_mask,
        )
        if self.encoding_ty == EncodingType.ONE_HOT:
            feat = torch.searchsorted(self.dist_bins.to(atom_ref_pos.device), feat)
        feat.masked_fill_(~mask, self.mask_value)

        return self.make_feature(feat.unsqueeze(-1))


@typecheck
def get_blocked_atom_pair_dists(
    positions: Float[Tensor, "b a 3"],
    atom_ref_space_uid: Int[Tensor, "b a"],
    q_idx: Int[Tensor, "bl bl_q"],
    kv_idx: Int[Tensor, "bl bl_kv"],
    block_atom_pair_mask: Bool[Tensor, "b bl bl_q bl_kv"],
) -> tuple[Float[Tensor, "b bl bl_q bl_kv"], Bool[Tensor, "b bl bl_q bl_kv"]]:
    q_pos = positions[:, q_idx]
    kv_pos = positions[:, kv_idx]

    blocked_pair_dists = cdist(q_pos, kv_pos)  # b bl bl_q bl_kv

    atom_ref_space_q = atom_ref_space_uid[:, q_idx]
    atom_ref_space_kv = atom_ref_space_uid[:, kv_idx]
    block_same_atom_ref_space = rearrange(
        atom_ref_space_q, "b bl a_q -> b bl a_q 1"
    ) == rearrange(atom_ref_space_kv, "b bl a_kv -> b bl 1 a_kv")

    block_atom_pair_mask &= block_same_atom_ref_space

    return blocked_pair_dists, block_atom_pair_mask

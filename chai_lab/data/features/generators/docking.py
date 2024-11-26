# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
import random
from dataclasses import dataclass
from typing import Any

import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.features.generators.token_pair_distance import TokenCenterDistance
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.model.utils import get_asym_id_from_subchain_id
from chai_lab.utils.defaults import default
from chai_lab.utils.tensor_utils import cdist, und, und_self
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class RestraintGroup:
    """
    Container for a docking constraint group --
    collection of chains with inter/intra distance constraints

    This class can be used to specify a set of chains to be
    grouped together for the docking feature
    """

    subchain_ids: list[str]
    noise_sigma: float
    dropout_prob: float
    atom_center_mask: list[Bool[Tensor, "_"]]
    atom_center_coords: list[Float[Tensor, "_ 3"]]

    def __post_init__(self) -> None:
        """Ensure params are consistent"""
        assert len(self.subchain_ids) == len(
            self.atom_center_coords
        ), f"{len(self.subchain_ids)=}, {len(self.atom_center_coords)=}"
        assert len(self.subchain_ids) == len(
            self.atom_center_mask
        ), f"{len(self.subchain_ids)=}, {len(self.atom_center_mask)=}"
        assert all(
            [
                len(mask) == len(coord)
                for coord, mask in zip(self.atom_center_coords, self.atom_center_mask)
            ]
        ), (
            f"{[len(x) for x in self.atom_center_coords]=}, "
            f"{[len(x) for x in self.atom_center_mask]=}"
        )

    def get_asym_ids(
        self,
        token_subchain_id: UInt8[Tensor, "n 4"],
        token_asym_id: Int[Tensor, "n"],
    ) -> list[int]:
        return [
            get_asym_id_from_subchain_id(
                subchain_id=subchain_id,
                source_pdb_chain_id=token_subchain_id,
                token_asym_id=token_asym_id,
            )
            for subchain_id in self.subchain_ids
        ]

    def __str__(self):
        return (
            f"{self.__class__.__name__}(subchain_ids={self.subchain_ids}, "
            f"atom_center_coords.shape={[x.shape for x in self.atom_center_coords]}, "
            f"atom_center_mask.shape={[x.shape for x in self.atom_center_mask]})"
        )


class DockingRestraintGenerator(FeatureGenerator):
    """Docking Feature Generator

    Works as follows:
        separate input chains into two groups by randomly
        partitioning asym_id's.
        Provide all token-center distances for chains within the
        same asm_id group.
        Mask all token-center distances for chains within the
        different asm_id groups.

    """

    def __init__(
        self,
        dist_bins: list[float] | None = None,
        coord_noise: tuple[float, float] = (0.0, 3.0),
        include_probability: float = 0.1,
        structure_dropout_prob: float = 0.0,
        chain_dropout_prob: float = 0.0,
        entity_types: list[EntityType] | None = None,
    ):
        dist_bins = dist_bins if dist_bins is not None else [0.0, 4.0, 8.0, 16.0]
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            encoding_ty=EncodingType.ONE_HOT,
            # one of dist_bins of rbf_radii is not None.
            num_classes=len(dist_bins) + 1,
            mult=1,
            can_mask=True,
        )
        self.token_dist_gen = TokenCenterDistance(dist_bins=dist_bins)

        # maintain consistent orders
        self.coord_noise = coord_noise
        self.include_probability = include_probability
        self.structure_dropout_prob = structure_dropout_prob
        self.chain_dropout_prob = chain_dropout_prob
        self.entity_types = set(
            [x.value for x in default(entity_types, [e for e in EntityType])]
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        maybe_constraint_dicts = batch["inputs"].get("docking_constraints", [[None]])[0]
        docking_constraints = batch["inputs"]["docking_constraints"] = (
            [RestraintGroup(**d) for d in maybe_constraint_dicts]
            if isinstance(maybe_constraint_dicts[0], dict)
            else None
        )

        return dict(
            all_atom_positions=batch["inputs"]["atom_gt_coords"],
            all_atom_mask=batch["inputs"]["atom_exists_mask"],
            token_single_mask=batch["inputs"]["token_exists_mask"],
            token_center_atom_index=batch["inputs"]["token_centre_atom_index"].long(),
            token_asym_id=batch["inputs"]["token_asym_id"].long(),
            token_subchain_id=batch["inputs"]["subchain_id"],
            token_entity_type=batch["inputs"]["token_entity_type"].long(),
            constraints=docking_constraints,
        )

    def apply_structure_dropout(
        self, feature: Tensor, prob: float | None = None
    ) -> Tensor:
        prob = default(prob, torch.rand(1).item())
        dropout_single_mask = torch.rand_like(feature.data[..., 0, 0].float()) < prob
        dropout_pair_mask = und_self(dropout_single_mask, "b i, b j -> b i j")
        feature = feature.masked_fill(dropout_pair_mask.unsqueeze(-1), self.mask_value)
        return feature

    def apply_chain_dropout(
        self, feature: Tensor, token_asym_id: Int[Tensor, "b n"]
    ) -> Tensor:
        structure_masks = []
        for i in range(token_asym_id.shape[0]):
            data_i, asym_i = feature.data[i], token_asym_id[i]
            unique_asyms = torch.unique(asym_i[asym_i != 0]).tolist()
            random.shuffle(unique_asyms)  # select chains to mask at random
            selected_asyms = unique_asyms[: random.randint(0, len(unique_asyms))]
            if len(selected_asyms) == 0:
                structure_masks.append(
                    torch.zeros_like(data_i[..., 0], dtype=torch.bool)
                )
                continue
            asyms_to_mask = torch.tensor(selected_asyms, device=data_i.device)
            asym_mask = torch.any(asym_i.unsqueeze(-1) == asyms_to_mask, dim=-1)
            structure_mask = und_self(asym_mask, "i, j -> i j")
            structure_masks.append(structure_mask)
        feature_mask = torch.stack(structure_masks, dim=0)
        feature = feature.masked_fill(feature_mask.unsqueeze(-1), self.mask_value)
        return feature

    @typecheck
    def _generate(
        self,
        all_atom_positions: Float[Tensor, "b a 3"],
        all_atom_mask: Bool[Tensor, "b a"],
        token_single_mask: Bool[Tensor, "b n"],
        token_center_atom_index: Int[Tensor, "b n"],
        token_asym_id: Int[Tensor, "b n"],
        token_entity_type: Int[Tensor, "b n"],
        token_subchain_id: UInt8[Tensor, "b n 4"],
        constraints: list[RestraintGroup] | None = None,
    ) -> Tensor:
        try:
            if constraints is not None:
                assert all_atom_positions.shape[0] == 1
                return self._generate_from_restraints(
                    token_asym_id=token_asym_id,
                    token_subchain_id=token_subchain_id,
                    constraints=constraints,
                )
        except Exception as e:
            logger.error(f"Error {e} generating docking constraints: {constraints}")

        return self._generate_from_batch(
            all_atom_positions=all_atom_positions,
            all_atom_mask=all_atom_mask,
            token_single_mask=token_single_mask,
            token_center_atom_index=token_center_atom_index,
            token_asym_id=token_asym_id,
            token_entity_type=token_entity_type,
        )

    def _asym_to_entity_type(
        self, asym_id: Int[Tensor, "n"], entity_type: Int[Tensor, "n"]
    ) -> dict[int, int]:
        unique_asyms: Tensor = torch.unique(asym_id[asym_id != 0])
        mapping = dict()
        for asym in unique_asyms.tolist():
            asym_mask = asym_id == asym
            mapping[int(asym)] = int(entity_type[asym_mask][0].item())
        return mapping

    @typecheck
    def _generate_from_batch(
        self,
        all_atom_positions=Float[Tensor, "b a 3"],
        all_atom_mask=Bool[Tensor, "b a"],
        token_single_mask=Bool[Tensor, "b n"],
        token_center_atom_index=Int[Tensor, "b n"],
        token_entity_type=Int[Tensor, "b n"],
        token_asym_id=Int[Tensor, "b n"],
    ) -> Tensor:
        sampled_noise = random.uniform(self.coord_noise[0], self.coord_noise[1])
        token_center_dists = self.token_dist_gen._generate(
            all_atom_positions=all_atom_positions
            + torch.randn_like(all_atom_positions) * sampled_noise,
            all_atom_mask=all_atom_mask,
            token_single_mask=token_single_mask,
            token_center_atom_index=token_center_atom_index,
        ).data
        for i in range(token_center_dists.shape[0]):
            asym_to_entity = self._asym_to_entity_type(
                token_asym_id[i], token_entity_type[i]
            )
            asym_include_list = [
                asym for asym, ety in asym_to_entity.items() if ety in self.entity_types
            ]
            asym_exclude_list = [
                asym
                for asym, ety in asym_to_entity.items()
                if ety not in self.entity_types
            ]
            # exclude other entity types
            asym_exclude_mask = torch.any(
                (token_asym_id[i].unsqueeze(-1) == torch.tensor(asym_exclude_list)),
                dim=-1,
            )
            token_center_dists[i, asym_exclude_mask] = self.mask_value
            token_center_dists[i, :, asym_exclude_mask] = self.mask_value
            if (
                random.random() < self.include_probability
                and len(asym_include_list) > 1
            ):
                # include distances between select chains
                random.shuffle(asym_include_list)
                partition_idx = random.randint(1, len(asym_include_list) - 1)
                _group_1, _group_2 = (
                    asym_include_list[:partition_idx],
                    asym_include_list[partition_idx:],
                )
                group_1, group_2 = torch.tensor(_group_1), torch.tensor(_group_2)
                # find positions of elements in first and second group
                group1_mask, group2_mask = [
                    torch.any((token_asym_id[i].unsqueeze(-1) == x), dim=-1)
                    for x in (group_1, group_2)
                ]
                partition_mask = und(group1_mask, group2_mask, "i, j -> i j")
                token_center_dists[i] = token_center_dists[i].masked_fill(
                    (partition_mask | partition_mask.T).unsqueeze(-1), self.mask_value
                )
            else:
                mask = torch.ones_like(token_center_dists[i], dtype=torch.bool)
                token_center_dists[i] = token_center_dists[i].masked_fill(
                    mask, self.mask_value
                )

        feature = self.make_feature(token_center_dists)
        if random.random() < self.structure_dropout_prob:
            feature = self.apply_structure_dropout(feature)
        elif random.random() < self.chain_dropout_prob:
            feature = self.apply_chain_dropout(feature, token_asym_id)
        return feature

    @typecheck
    def _generate_from_restraints(
        self,
        # constraints only supported with batch size 1
        token_asym_id: Int[Tensor, "1 n"],
        token_subchain_id: UInt8[Tensor, "1 n 4"],
        constraints: list[RestraintGroup],
    ) -> Tensor:
        logger.info(f"Generating docking feature from constraints: {constraints}")
        n, device = token_asym_id.shape[1], token_asym_id.device
        constraint_mat = torch.zeros(n, n, device=device, dtype=torch.float32)
        constraint_mask = torch.zeros(n, n, device=device, dtype=torch.bool)
        for constraint_group in constraints:
            # add constraints between members of each group
            coords = [
                x + torch.randn_like(x) * constraint_group.noise_sigma
                for x in constraint_group.atom_center_coords
            ]
            n_chains = len(constraint_group.subchain_ids)
            l_idx, r_idx = torch.triu_indices(n_chains, n_chains)
            chain_asyms = constraint_group.get_asym_ids(
                token_subchain_id=rearrange(token_subchain_id, "1 ... -> ..."),
                token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
            )
            for i, j in zip(l_idx.tolist(), r_idx.tolist()):
                constraint_mat, constraint_mask = self.add_restraint(
                    constraint_mat=constraint_mat,
                    constraint_mask=constraint_mask,
                    token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
                    chain1_asym_id=chain_asyms[i],
                    chain2_asym_id=chain_asyms[j],
                    chain1_coords=coords[i],
                    chain1_mask=constraint_group.atom_center_mask[i],
                    chain2_coords=coords[j],
                    chain2_mask=constraint_group.atom_center_mask[j],
                )
        # encode and apply mask
        feat = torch.searchsorted(
            self.token_dist_gen.dist_bins.to(constraint_mat.device), constraint_mat
        )
        feat = feat.masked_fill(~constraint_mask, self.mask_value)
        # add back batch dim
        constraint_mat = repeat(feat, "i j -> 1 i j 1")
        # apply structure dropout
        dropout = constraints[0].dropout_prob if len(constraints) > 0 else 0.0
        feature = self.make_feature(constraint_mat)
        feature = self.apply_structure_dropout(feature, prob=dropout)
        return feature

    @typecheck
    def add_restraint(
        self,
        constraint_mat: Float[Tensor, "n n"],
        constraint_mask: Bool[Tensor, "n n"],
        token_asym_id: Int[Tensor, "n"],
        chain1_asym_id: int,
        chain2_asym_id: int,
        chain1_coords: Float[Tensor, "c1 3"],
        chain2_coords: Float[Tensor, "c2 3"],
        chain1_mask: Bool[Tensor, "c1"],
        chain2_mask: Bool[Tensor, "c2"],
    ) -> tuple[Float[Tensor, "n n"], Bool[Tensor, "n n"]]:
        (c1_posns,) = torch.where(token_asym_id == chain1_asym_id)
        (c2_posns,) = torch.where(token_asym_id == chain2_asym_id)
        # make sure we have a coordinate for each position
        assert len(c1_posns) == len(
            chain1_coords
        ), f"{c1_posns.shape=}, {chain1_coords.shape=}"
        assert len(c2_posns) == len(
            chain2_coords
        ), f"{c2_posns.shape=}, {chain2_coords.shape=}"

        pairwise_dists = cdist(chain1_coords, chain2_coords)
        pairwise_mask = und(chain1_mask, chain2_mask, "i, j -> i j")
        pairwise_dists[~pairwise_mask] = -1.0
        # mask and fill the constraint matrix
        row_idxs = repeat(c1_posns, "i -> i c", c=len(c2_posns))
        col_idxs = repeat(c2_posns, "j -> r j", r=len(c1_posns))
        # fill constraints and mask
        constraint_mat[row_idxs, col_idxs] = pairwise_dists
        constraint_mat[col_idxs, row_idxs] = pairwise_dists
        constraint_mask[row_idxs, col_idxs] = pairwise_mask
        constraint_mask[col_idxs, row_idxs] = pairwise_mask
        return constraint_mat, constraint_mask

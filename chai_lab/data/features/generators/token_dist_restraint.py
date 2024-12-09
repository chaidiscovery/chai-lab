# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import dataclass

import numpy as np
import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.model.utils import get_asym_id_from_subchain_id
from chai_lab.utils.tensor_utils import tensorcode_to_string
from chai_lab.utils.typing import Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class RestraintGroup:
    """
    Container for a token pair distance restraint (contact)
    """

    left_residue_subchain_id: str
    right_residue_subchain_id: str
    left_residue_index: int
    right_residue_index: int
    right_residue_name: str
    left_residue_name: str
    distance_threshold: float

    def get_left_and_right_asym_ids(
        self,
        token_subchain_id: UInt8[Tensor, "n 4"],
        token_asym_id: Int[Tensor, "n"],
    ):
        left_asym_id = get_asym_id_from_subchain_id(
            subchain_id=self.left_residue_subchain_id,
            source_pdb_chain_id=token_subchain_id,
            token_asym_id=token_asym_id,
        )
        right_asym_id = get_asym_id_from_subchain_id(
            subchain_id=self.right_residue_subchain_id,
            source_pdb_chain_id=token_subchain_id,
            token_asym_id=token_asym_id,
        )
        return left_asym_id, right_asym_id

    def __str__(self):
        return (
            f"ConstraintGroup(left_residue_subchain_id={self.left_residue_subchain_id}, "
            f"right_residue_subchain_id={self.right_residue_subchain_id}, "
            f"left_residue_index={self.left_residue_index}, "
            f"right_residue_index={self.right_residue_index}, "
            f"right_residue_name={self.right_residue_name}, "
            f"left_residue_name={self.left_residue_name}, "
            f"distance_threshold={self.distance_threshold})"
        )


class TokenDistanceRestraint(FeatureGenerator):
    def __init__(
        self,
        include_probability: float = 1.0,
        size: int | float = 0.33,
        min_dist: int | float = 10.0,
        max_dist: int | float = 30.0,
        coord_noise: float = 0.0,
        num_rbf_radii: int = 5,
        query_entity_types: list[EntityType] | None = None,
        key_entity_types: list[EntityType] | None = None,
    ):
        """Randomly sample inter-chain token distance restraints

        Parameters:
            include_probability: Probability with which to include restraints
                for a given example. i.e. if include probability is 0.75, then 25%
                of the time, we will not sample any restraints for an example.
            size:  Number of restraints to sample. If 0 < size < 1, then the number
                of restraints will be determined as geom(size), independently for each
                example.
            min_dist: Minimum distance to encode restraints for
            max_dist: Maximum distance to encode restraints for
            coord_noise: gaussian noise with mean 0 and variance coord_noise
                added to coordinates before sampling restraints.
            num_rbf_radii: Number of radii to use for the radial basis function
                embedding of restraints
            query_entity_types: Entity types to consider when sampling "query" tokens
                for restraints. Defaults to all entity types.
            key_entity_types: Entity types to consider when sampling "key" tokens
                for restraints. Defaults to all entity types.

        NOTE: We only sample restraints between tokens if one of the tokens is in
            the query entity types and the other is in the key entity types.
        """
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            can_mask=False,
            encoding_ty=EncodingType.RBF,
            num_classes=num_rbf_radii,
            mult=1,
            ignore_index=-1.0,
        )
        self.ignore_idx = -1.0
        self.min_dist, self.max_dist = min_dist, max_dist
        self.coord_noise = coord_noise
        self.include_prob = include_probability
        self.size = size
        self.query_entity_types = torch.tensor(
            (
                [e.value for e in query_entity_types]
                if query_entity_types is not None
                else [e.value for e in EntityType]
            )
        ).long()
        self.key_entity_types = torch.tensor(
            [
                [e.value for e in key_entity_types]
                if key_entity_types is not None
                else [e.value for e in EntityType]
            ]
        ).long()

    def get_num_restraints(self, batch_size) -> list[int]:
        if 0 < self.size < 1:
            seles = np.random.geometric(self.size, size=batch_size)
            include_mask = np.random.uniform(size=batch_size) < self.include_prob
            seles[~include_mask] = 0
            return [int(x) for x in seles]
        return [int(self.size)] * batch_size

    def get_input_kwargs_from_batch(self, batch) -> dict:
        maybe_constraint_dicts = batch["inputs"].get("contact_constraints", [[None]])[0]
        contact_constraints = (
            [RestraintGroup(**d) for d in maybe_constraint_dicts]
            if isinstance(maybe_constraint_dicts[0], dict)
            else None
        )
        return dict(
            atom_gt_coords=batch["inputs"]["atom_gt_coords"],
            token_asym_id=batch["inputs"]["token_asym_id"].long(),
            token_residue_index=batch["inputs"]["token_residue_index"].long(),
            token_residue_names=batch["inputs"]["token_residue_name"],
            token_subchain_id=batch["inputs"]["subchain_id"],
            constraints=contact_constraints,
        )

    @typecheck
    def _generate(
        self,
        atom_gt_coords: Float[Tensor, "b a 3"],
        token_asym_id: Int[Tensor, "b n"],
        token_residue_index: Int[Tensor, "b n"],
        token_residue_names: UInt8[Tensor, "b n 8"],
        token_subchain_id: UInt8[Tensor, "b n 4"],
        constraints: list[RestraintGroup] | None = None,
    ) -> Float[Tensor, "b n n 1"]:
        try:
            if constraints is not None:
                assert atom_gt_coords.shape[0] == 1
                return self.generate_from_restraint(
                    token_asym_id=token_asym_id,
                    token_residue_index=token_residue_index,
                    token_residue_names=token_residue_names,
                    token_subchain_id=token_subchain_id,
                    constraints=constraints,
                )
        except Exception as e:
            logger.error(f"Error {e} generating distance constraints: {constraints}")

        batch, n = token_asym_id.shape
        device = token_asym_id.device
        constraint_mat = torch.full(
            (batch, n, n, 1),
            fill_value=self.ignore_idx,
            device=device,
            dtype=torch.float32,
        )
        return constraint_mat

    @typecheck
    def generate_from_restraint(
        self,
        token_asym_id: Int[Tensor, "1 n"],
        token_residue_index: Int[Tensor, "1 n"],
        token_residue_names: UInt8[Tensor, "1 n 8"],
        token_subchain_id: UInt8[Tensor, "1 n 4"],
        constraints: list[RestraintGroup],
    ) -> Tensor:
        logger.info(f"Generating distance feature from constraints: {constraints}")
        n, device = token_asym_id.shape[1], token_asym_id.device
        constraint_mat = torch.full(
            (n, n), fill_value=self.ignore_idx, device=device, dtype=torch.float32
        )
        for constraint_group in constraints:
            left_residue_asym_id, right_residue_asym_id = (
                constraint_group.get_left_and_right_asym_ids(
                    token_subchain_id=rearrange(token_subchain_id, "1 ... -> ..."),
                    token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
                )
            )
            constraint_mat = self.add_distance_restraint(
                constraint_mat=constraint_mat,
                token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
                token_residue_index=rearrange(token_residue_index, "1 ... -> ..."),
                token_residue_names=rearrange(token_residue_names, "1 ... -> ..."),
                left_residue_asym_id=left_residue_asym_id,
                right_residue_asym_id=right_residue_asym_id,
                left_residue_index=constraint_group.left_residue_index,
                right_residue_index=constraint_group.right_residue_index,
                right_residue_name=constraint_group.right_residue_name,
                left_residue_name=constraint_group.left_residue_name,
                distance_threshold=constraint_group.distance_threshold,
            )
        # encode and apply mask
        constraint_mat = repeat(constraint_mat, "i j -> 1 i j 1")
        return self.make_feature(constraint_mat)

    @typecheck
    def add_distance_restraint(
        self,
        constraint_mat: Float[Tensor, "n n"],
        token_asym_id: Int[Tensor, "n"],
        token_residue_index: Int[Tensor, "n"],
        token_residue_names: UInt8[Tensor, "n 8"],
        # asym id of the chain that binds in the pocket
        left_residue_asym_id: int,
        right_residue_asym_id: int,
        left_residue_index: int,
        right_residue_index: int,
        right_residue_name: str,
        left_residue_name: str,
        distance_threshold: float,
    ):
        left_asym_mask = token_asym_id == left_residue_asym_id
        right_asym_mask = token_asym_id == right_residue_asym_id
        left_index_mask = token_residue_index == left_residue_index
        right_index_mask = token_residue_index == right_residue_index
        left_residue_mask = left_asym_mask & left_index_mask
        right_residue_mask = right_asym_mask & right_index_mask
        # restraint should point to single residue pair
        assert torch.sum(left_residue_mask) == 1, (
            f"Expected unique residue but found {torch.sum(left_residue_mask)}\n"
            f"{left_residue_asym_id=}, {left_residue_index=}, {left_residue_name=}"
        )
        assert torch.sum(right_residue_mask) == 1, (
            f"Expected unique residue but found {torch.sum(right_residue_mask)}\n"
            f"{right_residue_asym_id=}, {right_residue_index=}, {right_residue_name=}"
        )
        # make sure the residue names in the constraint match the
        # ones we parsed
        left_res_name = token_residue_names[left_residue_mask]
        right_res_name = token_residue_names[right_residue_mask]
        expected_res_name = tensorcode_to_string(rearrange(left_res_name, "1 l -> l"))
        assert expected_res_name == left_residue_name, (
            f"Expected residue name {expected_res_name} but got " f"{left_residue_name}"
        )
        expected_res_name = tensorcode_to_string(rearrange(right_res_name, "1 l -> l"))
        assert expected_res_name == right_residue_name, (
            f"Expected residue name {expected_res_name} but got "
            f"{right_residue_name}"
        )
        # add constraint
        # NOTE: feature is *not* symmetric
        constraint_mat[left_residue_mask, right_residue_mask] = distance_threshold
        return constraint_mat

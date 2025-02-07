# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import dataclass

import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.features.generators.token_dist_restraint import (
    TokenDistanceRestraint,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.model.utils import get_asym_id_from_subchain_id
from chai_lab.utils.tensor_utils import tensorcode_to_string
from chai_lab.utils.typing import Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class RestraintGroup:
    """
    Container for a token pocket pair restraint group
    """

    # subchain ID of the pocket chain
    pocket_chain_subchain_id: str
    # subchain ID of the pocket token
    pocket_token_subchain_id: str
    # residue index of the pocket token
    pocket_token_residue_index: int
    # residue name of the pocket token
    pocket_token_residue_name: str
    # pocket distance threshold
    pocket_distance_threshold: float
    # optional subchain IDs

    def get_chain_and_token_asym_ids(
        self,
        token_subchain_id: UInt8[Tensor, "n 4"],
        token_asym_id: Int[Tensor, "n"],
    ):
        pocket_chain_asym_id = get_asym_id_from_subchain_id(
            subchain_id=self.pocket_chain_subchain_id,
            source_pdb_chain_id=token_subchain_id,
            token_asym_id=token_asym_id,
        )
        pocket_token_asym_id = get_asym_id_from_subchain_id(
            subchain_id=self.pocket_token_subchain_id,
            source_pdb_chain_id=token_subchain_id,
            token_asym_id=token_asym_id,
        )
        return pocket_chain_asym_id, pocket_token_asym_id

    def __str__(self):
        return (
            f"ConstraintGroup(pocket_chain_subchain_id={self.pocket_chain_subchain_id}, "
            f"pocket_token_subchain_id={self.pocket_token_subchain_id}, "
            f"pocket_token_residue_index={self.pocket_token_residue_index}, "
            f"pocket_token_residue_name={self.pocket_token_residue_name}, "
            f"pocket_distance_threshold={self.pocket_distance_threshold})"
        )


class TokenPairPocketRestraint(FeatureGenerator):
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
        """
        Derives pocket constraints by first generating pairwise distance restraints,
        and then selecting the query tokens that were assigned to some non-zero
        constraint.

        NOTE: Pocket restraints will only be sampled for tokens that are in the
            query entity types.
        """
        super().__init__(
            ty=FeatureType.TOKEN_PAIR,
            can_mask=False,
            encoding_ty=EncodingType.RBF,
            num_classes=num_rbf_radii,
            mult=1,
            ignore_index=-1.0,
        )
        # use distance restraint generator to sample pocket tokens/chains
        self.distance_restraint_gen = TokenDistanceRestraint(
            include_probability=include_probability,
            size=size,
            min_dist=min_dist,
            max_dist=max_dist,
            coord_noise=coord_noise,
            num_rbf_radii=num_rbf_radii,
            query_entity_types=query_entity_types,
            key_entity_types=key_entity_types,
        )
        self.ignore_idx = -1.0
        self.min_dist, self.max_dist = min_dist, max_dist
        self.coord_noise = coord_noise
        self.include_prob = include_probability
        self.size = size
        # override feature type
        self.ty = FeatureType.TOKEN_PAIR

    def get_input_kwargs_from_batch(self, batch) -> dict:
        # cast pocket constraints from dict back to dataclass
        maybe_constraint_dicts = batch["inputs"].get("pocket_constraints", [[None]])[0]
        pocket_constraints = batch["inputs"]["pocket_constraints"] = (
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
            constraints=pocket_constraints,
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
                return self.generate_from_restraints(
                    token_asym_id=token_asym_id,
                    token_residue_index=token_residue_index,
                    token_residue_names=token_residue_names,
                    token_subchain_id=token_subchain_id,
                    restraints=constraints,
                )
        except Exception as e:
            logger.error(f"Error {e} generating pocket constraints: {constraints}")

        # Return a null constraint matrix if constraints are not specified
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
    def generate_from_restraints(
        self,
        # only batch size 1 is supported
        token_asym_id: Int[Tensor, "1 n"],
        token_subchain_id: UInt8[Tensor, "1 n 4"],
        token_residue_index: Int[Tensor, "1 n"],
        token_residue_names: UInt8[Tensor, "1 n 8"],
        restraints: list[RestraintGroup],
    ) -> Tensor:
        logger.info(f"Generating pocket feature from constraints: {restraints}")
        n, device = token_asym_id.shape[1], token_asym_id.device
        constraint_mat = torch.full(
            (n, n), fill_value=self.ignore_idx, device=device, dtype=torch.float32
        )
        for restraint_group in restraints:
            pocket_chain_asym_id, pocket_token_asym_id = (
                restraint_group.get_chain_and_token_asym_ids(
                    token_subchain_id=rearrange(token_subchain_id, "1 ... -> ..."),
                    token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
                )
            )
            constraint_mat = self.add_pocket_restraint(
                restraint_mat=constraint_mat,
                token_asym_id=rearrange(token_asym_id, "1 ... -> ..."),
                token_residue_index=rearrange(token_residue_index, "1 ... -> ..."),
                token_residue_names=rearrange(token_residue_names, "1 ... -> ..."),
                pocket_chain_asym_id=pocket_chain_asym_id,
                pocket_token_asym_id=pocket_token_asym_id,
                pocket_token_residue_index=restraint_group.pocket_token_residue_index,
                pocket_token_residue_name=restraint_group.pocket_token_residue_name,
                pocket_distance_threshold=restraint_group.pocket_distance_threshold,
            )
        # encode and apply mask
        constraint_mat = repeat(constraint_mat, "i j -> 1 i j 1")
        return self.make_feature(constraint_mat)

    @typecheck
    def add_pocket_restraint(
        self,
        restraint_mat: Float[Tensor, "n n"],
        token_asym_id: Int[Tensor, "n"],
        token_residue_index: Int[Tensor, "n"],
        token_residue_names: UInt8[Tensor, "n 8"],
        # asym id of the chain that binds in the pocket
        pocket_chain_asym_id: int,
        # asym id of the token defining the pocket
        pocket_token_asym_id: int,
        # residue index of the pocket token
        pocket_token_residue_index: int,
        # residue name of the pocket token
        pocket_token_residue_name: str,
        # distance from the pocket token to pocket chain
        pocket_distance_threshold: float,
    ):
        pocket_chain_asym_mask = token_asym_id == pocket_chain_asym_id
        pocket_token_asym_mask = token_asym_id == pocket_token_asym_id
        pocket_token_residue_mask = token_residue_index == pocket_token_residue_index
        pocket_token_residue_mask &= pocket_token_asym_mask
        assert torch.sum(pocket_token_residue_mask) == 1, (
            f"Expected unique residue but found {torch.sum(pocket_token_residue_mask)}\n"
            f"{pocket_token_asym_id=}, {pocket_token_residue_index=}, "
            f"{pocket_token_residue_name=}"
        )
        pocket_token_res_name = token_residue_names[pocket_token_residue_mask]
        pocket_token_res_name = rearrange(pocket_token_res_name, "1 l -> l")
        expected_res_name = tensorcode_to_string(pocket_token_res_name)
        assert expected_res_name == pocket_token_residue_name, (
            f"Expected residue name {expected_res_name} but got "
            f"{pocket_token_residue_name}"
        )
        # add constraints between the pocket token and all other tokens in the pocket
        # chain
        # NOTE: feature is not symmetric
        restraint_mat[pocket_token_residue_mask, pocket_chain_asym_mask] = (
            pocket_distance_threshold
        )
        return restraint_mat

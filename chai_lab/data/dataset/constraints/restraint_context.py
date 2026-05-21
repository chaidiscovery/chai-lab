# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import asdict, dataclass

from torch import Tensor
from typing_extensions import assert_never

from chai_lab.data import residue_constants as rc
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.features.generators.docking import (
    RestraintGroup as DockingRestraint,
)
from chai_lab.data.features.generators.token_dist_restraint import (
    RestraintGroup as ContactRestraint,
)
from chai_lab.data.features.generators.token_pair_pocket_restraint import (
    RestraintGroup as PocketRestraint,
)
from chai_lab.data.parsing.restraints import (
    PairwiseInteraction,
    PairwiseInteractionType,
)
from chai_lab.utils.typing import typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class RestraintContext:
    docking_restraints: list[DockingRestraint] | None
    contact_restraints: list[ContactRestraint] | None
    pocket_restraints: list[PocketRestraint] | None

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"\n\tdocking_constraints {self.docking_restraints})"
            f"\n\tcontact_constraints {self.contact_restraints}"
            f"\n\tpocket_constraints {self.pocket_restraints}\n)"
        )

    def pad(self, *args, **kwargs) -> "RestraintContext":
        # No-op
        return RestraintContext(
            docking_restraints=self.docking_restraints,
            contact_restraints=self.contact_restraints,
            pocket_restraints=self.pocket_restraints,
        )

    def to_dict(self) -> dict[str, list[dict | None]]:
        return dict(
            docking_constraints=[asdict(c) for c in self.docking_restraints]
            if self.docking_restraints is not None
            else [None],
            contact_constraints=[asdict(c) for c in self.contact_restraints]
            if self.contact_restraints is not None
            else [None],
            pocket_constraints=[asdict(c) for c in self.pocket_restraints]
            if self.pocket_restraints is not None
            else [None],
        )

    @classmethod
    def empty(cls) -> "RestraintContext":
        return cls(
            docking_restraints=None,
            contact_restraints=None,
            pocket_restraints=None,
        )


def _is_cropped(
    chains: list[Chain],
    crop_idces: list[Tensor],
):
    return not all(
        [chain.num_tokens == len(crop) for chain, crop in zip(chains, crop_idces)]
    )


@typecheck
def load_manual_restraints_for_chai1(
    chains: list[Chain],
    crop_idces: list[Tensor] | None,
    provided_constraints: list[PairwiseInteraction],
) -> RestraintContext:
    """Load constraints from manual specification."""
    if len(chains) == 0 or len(provided_constraints) == 0:
        return RestraintContext.empty()
    assert crop_idces is None or not _is_cropped(chains, crop_idces)

    # For each of the constraints, add it into the constraint context
    docking_constraints: list[DockingRestraint] = []
    contact_constraints: list[ContactRestraint] = []
    pocket_constraints: list[PocketRestraint] = []

    logger.info(f"Loading {len(provided_constraints)} restraints...")
    for constraint in provided_constraints:
        match ctype := constraint.connection_type:
            case PairwiseInteractionType.COVALENT:
                # Covalent bonds are handled elsewhere, not as a constraint
                pass
            case PairwiseInteractionType.CONTACT:
                contact_parsed = ContactRestraint(
                    left_residue_subchain_id=constraint.chainA,
                    right_residue_subchain_id=constraint.chainB,
                    left_residue_index=constraint.res_idxA_pos - 1,  # 1 to 0-based
                    right_residue_index=constraint.res_idxB_pos - 1,  # 1 to 0-based
                    left_residue_name=rc.restype_1to3_with_x[constraint.res_idxA_name],
                    right_residue_name=rc.restype_1to3_with_x[constraint.res_idxB_name],
                    distance_threshold=constraint.max_dist_angstrom,
                )
                contact_constraints.append(contact_parsed)
            case PairwiseInteractionType.POCKET:
                # Treats A as the "chain level" and B as the "token level" constraints
                pocket_parsed = PocketRestraint(
                    pocket_chain_subchain_id=constraint.chainA,
                    pocket_token_subchain_id=constraint.chainB,
                    pocket_token_residue_index=constraint.res_idxB_pos - 1,
                    pocket_token_residue_name=rc.restype_1to3_with_x[
                        constraint.res_idxB_name
                    ],
                    pocket_distance_threshold=constraint.max_dist_angstrom,
                )
                pocket_constraints.append(pocket_parsed)
            case _:
                assert_never(ctype)
    return RestraintContext(
        docking_restraints=docking_constraints if docking_constraints else None,
        contact_restraints=contact_constraints if contact_constraints else None,
        pocket_restraints=pocket_constraints if pocket_constraints else None,
    )

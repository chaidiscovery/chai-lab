# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

import logging
from dataclasses import asdict, dataclass
from typing import Any, assert_never

from torch import Tensor

from chai_lab.data import residue_constants as rc
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.features.generators.docking import (
    ConstraintGroup as DockingConstraint,
)
from chai_lab.data.features.generators.token_dist_restraint import (
    ConstraintGroup as ContactConstraint,
)
from chai_lab.data.features.generators.token_pair_pocket_restraint import (
    ConstraintGroup as PocketConstraint,
)
from chai_lab.data.parsing.constraints import (
    PairwiseInteraction,
    PairwiseInteractionType,
)
from chai_lab.utils.typing import typecheck


@typecheck
@dataclass
class ConstraintContext:
    docking_constraints: list[DockingConstraint] | None
    contact_constraints: list[ContactConstraint] | None
    pocket_constraints: list[PocketConstraint] | None

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"\n\tdocking_constraints {self.docking_constraints})"
            f"\n\tcontact_constraints {self.contact_constraints}"
            f"\n\tpocket_constraints {self.pocket_constraints}\n)"
        )

    def pad(self, *args, **kwargs) -> "ConstraintContext":
        # No-op
        return ConstraintContext(
            docking_constraints=self.docking_constraints,
            contact_constraints=self.contact_constraints,
            pocket_constraints=self.pocket_constraints,
        )

    def to_dict(self) -> dict[str, list[dict | Any]]:
        return dict(
            docking_constraints=[asdict(c) for c in self.docking_constraints]
            if self.docking_constraints is not None
            else [None],
            contact_constraints=[asdict(c) for c in self.contact_constraints]
            if self.contact_constraints is not None
            else [None],
            pocket_constraints=[asdict(c) for c in self.pocket_constraints]
            if self.pocket_constraints is not None
            else [None],
        )

    @classmethod
    def empty(cls) -> "ConstraintContext":
        return cls(
            docking_constraints=None,
            contact_constraints=None,
            pocket_constraints=None,
        )


def _is_cropped(
    chains: list[Chain],
    crop_idces: list[Tensor],
):
    return not all(
        [chain.num_tokens == len(crop) for chain, crop in zip(chains, crop_idces)]
    )


@typecheck
def load_manual_constraints_for_chai1(
    chains: list[Chain],
    crop_idces: list[Tensor] | None,
    provided_constraints: list[PairwiseInteraction],
) -> ConstraintContext:
    """Load constraints from manual specification."""
    if len(chains) == 0 or len(provided_constraints) == 0:
        return ConstraintContext.empty()
    assert crop_idces is None or not _is_cropped(chains, crop_idces)

    # For each of the constraints, add it into the constraint context
    docking_constraints: list[DockingConstraint] = []
    contact_constraints: list[ContactConstraint] = []
    pocket_constraints: list[PocketConstraint] = []

    logging.info(f"Loading {len(provided_constraints)} constraints...")
    for constraint in provided_constraints:
        match ctype := constraint.connection_type:
            case PairwiseInteractionType.COVALENT:
                # Covalent bonds are handled elsewhere, not as a constraint
                pass
            case PairwiseInteractionType.CONTACT:
                contact_parsed = ContactConstraint(
                    left_residue_subchain_id=constraint.chainA,
                    right_residue_subchain_id=constraint.chainB,
                    left_residue_index=constraint.res_idxA_pos,
                    right_residue_index=constraint.res_idxB_pos,
                    left_residue_name=rc.restype_1to3_with_x[constraint.res_idxA_name],
                    right_residue_name=rc.restype_1to3_with_x[constraint.res_idxB_name],
                    distance_threshold=constraint.max_dist_angstrom,
                )
                contact_constraints.append(contact_parsed)
            case PairwiseInteractionType.POCKET:
                # Treats A as the "chain level" and B as the "token level" constraints
                pocket_parsed = PocketConstraint(
                    pocket_chain_subchain_id=constraint.chainA,
                    pocket_token_subchain_id=constraint.chainB,
                    pocket_token_residue_index=constraint.res_idxB_pos,
                    pocket_token_residue_name=rc.restype_1to3_with_x[
                        constraint.res_idxB_name
                    ],
                    pocket_distance_threshold=constraint.max_dist_angstrom,
                )
                pocket_constraints.append(pocket_parsed)
            case _:
                assert_never(ctype)
    return ConstraintContext(
        docking_constraints=docking_constraints if docking_constraints else None,
        contact_constraints=contact_constraints if contact_constraints else None,
        pocket_constraints=pocket_constraints if pocket_constraints else None,
    )

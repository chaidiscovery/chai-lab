# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

from dataclasses import dataclass

from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.parsing.structure.all_atom_entity_data import AllAtomEntityData


@dataclass
class Chain:
    # The untokenized entity data
    entity_data: AllAtomEntityData

    # The tokenized chain, derived from the entity data
    structure_context: AllAtomStructureContext

    def __str__(self) -> str:
        return f"{self.__class__.__name__}(entity_data={self.entity_data})"

    @property
    def num_tokens(self) -> int:
        return self.structure_context.num_tokens

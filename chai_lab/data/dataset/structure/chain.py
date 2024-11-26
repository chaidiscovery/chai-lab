# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

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

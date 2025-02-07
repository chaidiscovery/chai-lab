# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from enum import Enum

import gemmi

logger = logging.getLogger(__name__)


class EntityType(Enum):
    PROTEIN = 0
    RNA = 1
    DNA = 2
    LIGAND = 3
    POLYMER_HYBRID = 4
    WATER = 5
    UNKNOWN = 6
    MANUAL_GLYCAN = 7  # NOTE glycan parsing


def get_entity_type(entity: gemmi.Entity) -> EntityType:
    match (entity.entity_type, entity.polymer_type):
        case (
            gemmi.EntityType.Polymer,
            gemmi.PolymerType.PeptideL
            | gemmi.PolymerType.PeptideD,
        ):
            return EntityType.PROTEIN
        case (gemmi.EntityType.Polymer, gemmi.PolymerType.Dna):
            return EntityType.DNA
        case (gemmi.EntityType.Polymer, gemmi.PolymerType.Rna):
            return EntityType.RNA
        case (gemmi.EntityType.Polymer, gemmi.PolymerType.DnaRnaHybrid):
            return EntityType.POLYMER_HYBRID
        case (gemmi.EntityType.NonPolymer, _):
            return EntityType.LIGAND
        case (gemmi.EntityType.Branched, _):
            # multi residue ligand; do NOT return MANUAL_GLYCAN here
            return EntityType.LIGAND
        case (gemmi.EntityType.Water, _):
            return EntityType.WATER
        case (gemmi.EntityType.Branched, _):
            return EntityType.LIGAND
        case _:
            return EntityType.UNKNOWN

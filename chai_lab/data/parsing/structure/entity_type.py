# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from enum import Enum

logger = logging.getLogger(__name__)


class EntityType(Enum):
    PROTEIN = 0
    RNA = 1
    DNA = 2
    LIGAND = 3
    POLYMER_HYBRID = 4
    WATER = 5
    UNKNOWN = 6

# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

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

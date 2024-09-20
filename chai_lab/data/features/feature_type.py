# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

from enum import Enum


class FeatureType(Enum):
    RESIDUE = "RESIDUE"
    PAIR = "PAIR"
    MSA = "MSA"
    TEMPLATES = "TEMPLATES"
    TOKEN = "TOKEN"
    TOKEN_PAIR = "TOKEN_PAIR"
    ATOM = "ATOM"
    ATOM_PAIR = "ATOM_PAIR"

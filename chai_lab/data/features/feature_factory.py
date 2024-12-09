# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Helper methods for generating model input features"""

import logging

from torch import Tensor

from chai_lab.data.features.generators.base import FeatureGenerator

logger = logging.getLogger(__name__)


class FeatureFactory:
    generators: dict[str, FeatureGenerator]

    def __init__(self, generators: dict[str, FeatureGenerator]):
        self.generators = generators

    def generate(self, batch) -> dict[str, Tensor]:
        return {name: gen.generate(batch) for name, gen in self.generators.items()}

    def __repr__(self) -> str:
        return f"Feature factory, {len(self.generators)=}"

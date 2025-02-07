# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Feature Generator ABC and Default implementation"""

from abc import ABC
from enum import Enum

import torch
from beartype import beartype as typechecker
from torch import Tensor
from typing_extensions import assert_never

from chai_lab.data.features.feature_type import FeatureType


class EncodingType(Enum):
    ONE_HOT = "one-hot"
    RBF = "rbf"
    FOURIER = "fourier"
    IDENTITY = "identity"
    ESM = "esm"
    OUTERSUM = "outersum"


def cast_feature(
    feature: Tensor,
    encoding_ty: EncodingType,
):
    match encoding_ty:
        case EncodingType.IDENTITY:
            feature = feature.float()
            # safety check
            assert feature.abs().max() < 100, feature
            return feature
        case EncodingType.RBF | EncodingType.FOURIER:
            assert feature.dtype in (torch.float16, torch.float32, torch.bfloat16)
            return feature
        case EncodingType.ONE_HOT | EncodingType.OUTERSUM:
            if feature.dtype not in {
                torch.long,
                torch.int,
                torch.int16,
                torch.int8,
                torch.uint8,
            }:
                raise ValueError(
                    f"dtype {feature.dtype} is not a valid type for {encoding_ty}"
                )
            return feature
        case EncodingType.ESM:
            return feature

    assert_never(encoding_ty)  # Enum exhaustiveness check


class FeatureGenerator(ABC):
    @typechecker
    def __init__(
        self,
        ty: FeatureType,
        encoding_ty: EncodingType,
        num_classes: int = -1,
        mult: int = 1,
        ignore_index: float = -100.0,
        can_mask: bool = True,  # marks existing, but unknown values (e.g. atom position)
    ):
        self.ty = ty
        self.encoding_ty = encoding_ty
        self.num_classes = num_classes
        self.mult = mult
        self.ignore_index = ignore_index
        self.can_mask = can_mask

    @property
    def mask_value(self) -> int | float | Tensor:
        """Get value used to mask this feature"""
        match self.encoding_ty:
            case EncodingType.ONE_HOT | EncodingType.OUTERSUM:
                return self.num_classes
            case EncodingType.FOURIER | EncodingType.RBF:
                return -100.0
            case EncodingType.IDENTITY:
                assert self.can_mask
                mask = torch.zeros(self.num_classes + int(self.can_mask))
                mask[-1] = 1  # last channel is 1 for masked-out items
                return mask
            case EncodingType.ESM:
                return 0.0

        assert_never(self.encoding_ty)  # Enum exhaustiveness check

    def generate(self, batch) -> Tensor:
        """Generate a feature"""
        kwargs = self.get_input_kwargs_from_batch(batch)
        feature = self._generate(**kwargs)
        return feature

    def _generate(self, *args, **kwargs) -> Tensor:
        """Generate a feature"""
        raise NotImplementedError("implement me")

    def get_input_kwargs_from_batch(self, batch) -> dict:
        """Get input keyword arguments to pass to _generate"""
        raise NotImplementedError("implement me")

    def make_feature(self, data: Tensor) -> Tensor:
        """Checks and converts dtype if necessary"""
        return cast_feature(data, encoding_ty=self.encoding_ty)

    def __repr__(self):
        return f"[FeatureGenerator] : type: {self.ty}"

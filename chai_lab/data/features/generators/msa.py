# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import Any

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.dataset.msas.msa_context import NO_PAIRING_KEY
from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.parsing.msas.data_source import msa_dataset_source_to_int
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order
from chai_lab.utils.tensor_utils import masked_mean
from chai_lab.utils.typing import Bool, Int, UInt8, typecheck


class MSAFeatureGenerator(FeatureGenerator):
    """Generates feature for one-hot encoding of processed MSA, same classes as restype."""

    def __init__(self):
        num_res_ty = len(residue_types_with_nucleotides_order)
        super().__init__(
            ty=FeatureType.MSA,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=False,
            num_classes=num_res_ty,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            msa_tokens=batch["inputs"]["msa_tokens"],
        )

    @typecheck
    def _generate(
        self,
        msa_tokens: UInt8[Tensor, "batch depth tokens"],
    ) -> Tensor:
        """Generate based on an input of one-hot encoded MSA"""
        return self.make_feature(data=msa_tokens.unsqueeze(-1))


class MSAHasDeletionGenerator(FeatureGenerator):
    """Binary feature for if there is a deletion to the left of each position."""

    def __init__(self):
        super().__init__(
            ty=FeatureType.MSA,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(msa_deletion_matrix=batch["inputs"]["msa_deletion_matrix"])

    @typecheck
    def _generate(
        self,
        msa_deletion_matrix: UInt8[Tensor, "batch depth tokens"],
    ) -> Tensor:
        has_deletion = msa_deletion_matrix > 0
        return self.make_feature(data=has_deletion.unsqueeze(-1))


class MSADeletionValueGenerator(FeatureGenerator):
    """Raw deletion counts left of the current position, with addtional scaling.
    Scaling is given by s(d) = 2 / pi * arctan(d / 3)
    """

    def __init__(self):
        super().__init__(
            ty=FeatureType.MSA,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(msa_deletion_matrix=batch["inputs"]["msa_deletion_matrix"])

    @typecheck
    def _generate(
        self,
        msa_deletion_matrix: UInt8[Tensor, "batch depth tokens"],
    ) -> Tensor:
        d_scaled = 2.0 / torch.pi * torch.arctan(msa_deletion_matrix.float() / 3.0)
        return self.make_feature(data=d_scaled.unsqueeze(-1))


class MSAProfileGenerator(FeatureGenerator):
    """MSA profile - distribution across residue types BEFORE processing"""

    def __init__(self):
        self.num_res_ty = len(residue_types_with_nucleotides_order)
        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=self.num_res_ty,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            main_msa_tokens=batch["inputs"]["main_msa_tokens"],
            main_msa_mask=batch["inputs"]["main_msa_mask"],
        )

    @typecheck
    def _generate(
        self,
        main_msa_tokens: UInt8[Tensor, "batch depth tokens"],
        main_msa_mask: Bool[Tensor, "batch depth tokens"],
    ) -> Tensor:
        """Optimized implementation based on torch.scatter_add"""
        batch, _, tokens = main_msa_tokens.shape

        unnormalized_profile = torch.zeros(
            (batch, tokens, self.num_res_ty), dtype=main_msa_tokens.dtype
        ).scatter_add(
            dim=2,
            index=rearrange(
                main_msa_tokens.long(), "batch depth tokens -> batch tokens depth"
            ),
            src=rearrange(
                main_msa_mask.to(main_msa_tokens.dtype),
                "batch depth tokens -> batch tokens depth",
            ),
        )
        denom = unnormalized_profile.sum(dim=-1, keepdim=True).clamp_min_(1)
        profile = unnormalized_profile / denom

        return self.make_feature(data=profile)


class MSADeletionMeanGenerator(FeatureGenerator):
    """MSA deletion mean - mean number of deletions at each position in main MSA."""

    def __init__(self):
        super().__init__(
            ty=FeatureType.TOKEN,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            main_msa_mask=batch["inputs"]["main_msa_mask"],
            main_msa_deletion_matrix=batch["inputs"]["main_msa_deletion_matrix"],
        )

    @typecheck
    def _generate(
        self,
        main_msa_mask: Bool[Tensor, "batch depth tokens"],
        main_msa_deletion_matrix: UInt8[Tensor, "batch depth tokens"],
    ) -> Tensor:
        """Mean number of deletions at each position in main MSA."""
        # Average out the depth to get per-tokens
        mean_deletion_matrix = masked_mean(
            mask=main_msa_mask, value=main_msa_deletion_matrix.float(), dim=1
        )
        return self.make_feature(data=mean_deletion_matrix.unsqueeze(-1))


class IsPairedMSAGenerator(FeatureGenerator):
    """
    Assuming pairkey is species, for each token informs if it comes
    from the same species as the first token in complex.
    """

    def __init__(self):
        super().__init__(
            ty=FeatureType.MSA,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=1,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            msa_mask=batch["inputs"]["msa_mask"],
            pairing_key=batch["inputs"]["msa_pairkey"],
        )

    @typecheck
    def _generate(
        self,
        msa_mask: Bool[Tensor, "batch depth tokens"],
        pairing_key: Int[Tensor, "batch depth tokens"],
    ) -> Tensor:
        can_be_paired = msa_mask & (pairing_key != NO_PAIRING_KEY)

        is_paired = pairing_key == pairing_key[..., :1]
        is_paired = (is_paired & can_be_paired).to(torch.uint8)

        return self.make_feature(data=is_paired.unsqueeze(-1))


class MSADataSourceGenerator(FeatureGenerator):
    """
    MSA data source for each MSA token
    """

    def __init__(
        self,
        num_classes: int = 6,  # chai1 : 5 classes + mask val
    ):
        super().__init__(
            ty=FeatureType.MSA,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=num_classes,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            msa_mask=batch["inputs"]["msa_mask"],
            msa_sequence_source=batch["inputs"]["msa_sequence_source"],
        )

    @typecheck
    def _generate(
        self,
        msa_mask: Bool[Tensor, "batch depth tokens"],
        msa_sequence_source: UInt8[Tensor, "batch depth tokens"],
    ) -> Tensor:
        from chai_lab.data.parsing.msas.data_source import MSADataSource

        query = msa_dataset_source_to_int[MSADataSource.QUERY]
        none = msa_dataset_source_to_int[MSADataSource.NONE]
        # chai-1 specific: replace QUERY with NONE
        msa_sequence_source[msa_sequence_source.eq(query)] = none
        # use none for masking.
        msa_sequence_source = msa_sequence_source.masked_fill(~msa_mask, none)
        return self.make_feature(data=msa_sequence_source.unsqueeze(-1))

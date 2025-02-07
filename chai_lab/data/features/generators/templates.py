# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Feature generators for templates. This includes the following:
- Template mask (includes both the psuedo beta mask and backbone frame mask)
- Template unit vector generator
- Template residue type generator
- Template distogram generator
"""

import logging
from typing import Any

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


class TemplateMaskGenerator(FeatureGenerator):
    def __init__(self):
        super().__init__(
            ty=FeatureType.TEMPLATES,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=2,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            template_backbone_frame_mask=batch["inputs"][
                "template_backbone_frame_mask"
            ],
            template_pseudo_beta_mask=batch["inputs"]["template_pseudo_beta_mask"],
            asym_ids=batch["inputs"]["token_asym_id"].type(torch.int32),
        )

    def _generate(
        self,
        template_backbone_frame_mask: Bool[Tensor, "batch templ tokens"],
        template_pseudo_beta_mask: Bool[Tensor, "batch templ tokens"],
        asym_ids: Int[Tensor, "batch tokens"],
    ) -> Tensor:
        same_asym = rearrange(asym_ids, "b t -> b 1 t 1 1") == rearrange(
            asym_ids, "b t -> b 1 1 t 1"
        )
        # Line 1: backbone frame mask
        # (b t n n)
        bij_backbone = rearrange(
            template_backbone_frame_mask, "b t n -> b t n 1 1"
        ) * rearrange(template_backbone_frame_mask, "b t n -> b t 1 n 1")

        # Line 2: backbone pseudo beta mask
        # (b t n n)
        bij_pseudo_beta = rearrange(
            template_pseudo_beta_mask, "b t n -> b t n 1 1"
        ) * rearrange(template_pseudo_beta_mask, "b t n -> b t 1 n 1")

        mask_feat = torch.cat([bij_backbone, bij_pseudo_beta], dim=-1).float()

        return self.make_feature(mask_feat.float() * same_asym.float())


class TemplateUnitVectorGenerator(FeatureGenerator):
    """Generates feature for template unit vector"""

    def __init__(self):
        super().__init__(
            ty=FeatureType.TEMPLATES,
            encoding_ty=EncodingType.IDENTITY,
            can_mask=False,
            num_classes=3,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            template_unit_vector=batch["inputs"]["template_unit_vector"],
            asym_ids=batch["inputs"]["token_asym_id"].to(torch.int32),
        )

    @typecheck
    def _generate(
        self,
        template_unit_vector: Float[Tensor, "batch templ tokens tokens 3"],
        asym_ids: Int[Tensor, "batch tokens"],
    ) -> Tensor:
        same_asym = rearrange(asym_ids, "b t -> b 1 t 1 1") == rearrange(
            asym_ids, "b t -> b 1 1 t 1"
        )
        same_asym = same_asym.to(template_unit_vector.dtype)
        # mask out pairs with different asyms
        template_unit_vector = template_unit_vector * same_asym
        return self.make_feature(template_unit_vector)


class TemplateResTypeGenerator(FeatureGenerator):
    """Generates feature for one-hot encoding of templates, same classes as restype."""

    def __init__(self, embed_dim=32):
        num_res_ty = len(residue_types_with_nucleotides_order)
        super().__init__(
            ty=FeatureType.TEMPLATES,
            encoding_ty=EncodingType.OUTERSUM,
            can_mask=False,
            num_classes=num_res_ty,
            mult=1,
        )

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            template_tokens=batch["inputs"]["template_restype"].type(torch.uint8),
        )

    @typecheck
    def _generate(
        self,
        template_tokens: UInt8[Tensor, "batch templ tokens"],
    ) -> Tensor:
        return self.make_feature(data=template_tokens.unsqueeze(-1))


class TemplateDistogramGenerator(FeatureGenerator):
    """Generates feature for distogram of templates."""

    def __init__(
        self,
        min_dist_bin: float = 3.25,
        max_dist_bin: float = 50.75,
        n_dist_bin: int = 38,
    ):
        super().__init__(
            ty=FeatureType.TEMPLATES,
            encoding_ty=EncodingType.ONE_HOT,
            can_mask=True,
            num_classes=n_dist_bin,
            mult=1,
        )
        self.dist_bins = torch.linspace(min_dist_bin, max_dist_bin, n_dist_bin)[1:]

    def get_input_kwargs_from_batch(self, batch: dict[str, Any]) -> dict:
        return dict(
            template_distances=batch["inputs"]["template_distances"],
            asym_ids=batch["inputs"]["token_asym_id"].to(torch.int32),
        )

    @typecheck
    def _generate(
        self,
        template_distances: Float[Tensor, "batch templ tokens tokens"],
        asym_ids: Int[Tensor, "batch tokens"],
    ) -> Tensor:
        discretized = torch.searchsorted(self.dist_bins, template_distances)
        same_asym = rearrange(asym_ids, "b t -> b 1 t 1") == rearrange(
            asym_ids, "b t -> b 1 1 t"
        )
        discretized = torch.masked_fill(discretized, ~same_asym, self.mask_value)
        return self.make_feature(data=discretized.unsqueeze(-1))

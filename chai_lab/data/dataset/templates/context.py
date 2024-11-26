# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import asdict, dataclass

import torch
from torch import Tensor
from torch.nn import functional as F

from chai_lab.data import residue_constants as rc
from chai_lab.utils.defaults import default
from chai_lab.utils.typing import Bool, Float, Int, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass(frozen=True)
class TemplateContext:
    """Context for templates; always aligned by construction."""

    template_restype: Int[Tensor, "n_templates n_tokens"]
    template_pseudo_beta_mask: Bool[Tensor, "n_templates n_tokens"]
    template_backbone_frame_mask: Bool[Tensor, "n_templates n_tokens"]
    template_distances: Float[Tensor, "n_templates n_tokens n_tokens"]
    template_unit_vector: Float[Tensor, "n_templates n_tokens n_tokens 3"]

    def __str__(self) -> str:
        return (
            f"TemplateContext(num_templates={self.num_templates}, "
            f"num_nonnull_templates={self.num_nonnull_templates}, "
            f"num_tokens={self.num_tokens})"
        )

    @property
    def num_tokens(self) -> int:
        return self.template_restype.shape[1]

    @property
    def num_templates(self) -> int:
        return self.template_restype.shape[0]

    @property
    def num_nonnull_templates(self) -> int:
        """Number of templates that aren't all null padding templates."""
        template_exists = self.template_mask.any(dim=-1).int()
        return int(template_exists.sum().item())

    @property
    def template_mask(self) -> Bool[Tensor, "n_templates n_tokens"]:
        return self.template_restype != rc.residue_types_with_nucleotides_order["-"]

    def to_dict(self) -> dict[str, torch.Tensor]:
        retval = asdict(self)
        retval.update(
            {
                "num_templates": torch.tensor(self.num_nonnull_templates),
                "template_mask": self.template_mask,
            }
        )
        return retval

    @classmethod
    def empty(cls, n_templates: int, n_tokens: int) -> "TemplateContext":
        return cls(
            template_restype=torch.full(
                (n_templates, n_tokens),
                fill_value=rc.residue_types_with_nucleotides_order["-"],
                dtype=torch.int32,
            ),
            template_pseudo_beta_mask=torch.zeros(
                n_templates, n_tokens, dtype=torch.bool
            ),
            template_backbone_frame_mask=torch.zeros(
                n_templates, n_tokens, dtype=torch.bool
            ),
            template_distances=torch.zeros(
                n_templates, n_tokens, n_tokens, dtype=torch.float32
            ),
            template_unit_vector=torch.zeros(
                n_templates, n_tokens, n_tokens, 3, dtype=torch.float32
            ),
        )

    def index_select(self, idxs: Int[Tensor, "n"]) -> "TemplateContext":
        return TemplateContext(
            template_restype=self.template_restype[:, idxs],
            template_pseudo_beta_mask=self.template_pseudo_beta_mask[:, idxs],
            template_backbone_frame_mask=self.template_backbone_frame_mask[:, idxs],
            template_distances=self.template_distances[:, idxs][:, :, idxs],
            template_unit_vector=self.template_unit_vector[:, idxs][:, :, idxs],
        )

    # @classmethod
    # def merge(
    #     cls,
    #     templates: list["TemplateContext"],
    # ) -> "TemplateContext":
    #     """Merge template contexts along the template dimensions."""
    #     # n_token can be simply concatenated
    #     logger.debug(f"Merging {len(templates)} templates")

    #     # Handle case where we get an empty list (no templates to merge)
    #     if len(templates) == 0:
    #         return cls.empty(n_templates=4, n_tokens=1)

    #     # Pad each template_restype's template_dimension to match the largest
    #     # NOTE count num_templates here, NOT num_nonnull_templates
    #     n_templates_new: int = max(t.num_templates for t in templates)
    #     padded_templates = [t.pad(max_templates=n_templates_new) for t in templates]
    #     new_template_restype = torch.cat(
    #         [t.template_restype for t in padded_templates],
    #         dim=1,  # Concat on sequence dim
    #     )
    #     new_template_pseudo_beta_mask = torch.cat(
    #         [t.template_pseudo_beta_mask for t in padded_templates],
    #         dim=1,
    #     )
    #     new_template_backbone_frame_mask = torch.cat(
    #         [t.template_backbone_frame_mask for t in padded_templates],
    #         dim=1,
    #     )

    #     # Number of tokens after concatenation along token dim
    #     n_token_new = new_template_restype.shape[1]

    #     # n_token x n_token must be tiled into a square matrix
    #     # These indices like [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, ...] indicate the region
    #     # of the square matrix that corresponds to each template.
    #     template_indices = torch.repeat_interleave(
    #         input=torch.arange(len(templates), device=new_template_restype.device),
    #         repeats=torch.tensor([t.template_restype.shape[-1] for t in templates]),
    #     )
    #     # Sample template and token dim
    #     assert template_indices.shape[0] == n_token_new

    #     new_template_distances = torch.zeros(
    #         n_templates_new, n_token_new, n_token_new, dtype=torch.float32
    #     )
    #     new_template_unit_vector = torch.zeros(
    #         n_templates_new, n_token_new, n_token_new, 3, dtype=torch.float32
    #     )

    #     # For each template, find the block that it corresponds to and copy in the data
    #     for i, t in enumerate(templates):
    #         m = template_indices == i
    #         mask = m[:, None] * m[None, :]
    #         idx = torch.arange(t.template_distances.shape[0])
    #         new_template_distances[idx.unsqueeze(1), mask] = (
    #             t.template_distances.flatten(1, 2)
    #         )
    #         new_template_unit_vector[idx.unsqueeze(1), mask] = (
    #             t.template_unit_vector.flatten(1, 2)
    #         )

    #     return cls(
    #         template_restype=new_template_restype,
    #         template_pseudo_beta_mask=new_template_pseudo_beta_mask,
    #         template_backbone_frame_mask=new_template_backbone_frame_mask,
    #         template_distances=new_template_distances,
    #         template_unit_vector=new_template_unit_vector,
    #     )

    def pad(
        self,
        max_templates: int | None = None,
        max_tokens: int | None = None,
    ) -> "TemplateContext":
        """Pad to the given number of templates and tokens."""
        max_templates = default(max_templates, self.num_templates)
        assert (
            self.num_templates <= max_templates
        ), f"Cannot pad templates containing {self.num_templates} templates to {max_templates} templates"
        n_pad_templates = max_templates - self.num_templates

        max_tokens = default(max_tokens, self.num_tokens)
        assert (
            self.num_tokens <= max_tokens
        ), f"Cannot pad templates containing {self.num_tokens} tokens to {max_tokens} tokens"
        n_pad_tokens = max_tokens - self.num_tokens

        if n_pad_templates == 0 and n_pad_tokens == 0:  # Exact match yay
            return self

        logger.debug(f"Padding templates by {n_pad_templates=} {n_pad_tokens=}")

        # Padding works from last dim forward in pairs of padding (left, right)
        # - (0, n_pad_tokens) = pad nothing on left, pad by n_pad_tokens on right for
        #   last dim
        # - (0, 0, 0, n_pad_tokens, 0, n_pad_tokens) = pad nothing on last dim, but pad
        #   next two dims
        pad_dims_template = (0, n_pad_templates)
        pad_dims_token = (0, n_pad_tokens)
        return TemplateContext(
            template_restype=F.pad(
                self.template_restype,
                pad=pad_dims_token + pad_dims_template,
                value=rc.residue_types_with_nucleotides_order["-"],
            ),
            template_pseudo_beta_mask=F.pad(
                self.template_pseudo_beta_mask,
                pad=pad_dims_token + pad_dims_template,
            ),
            template_backbone_frame_mask=F.pad(
                self.template_backbone_frame_mask,
                pad=pad_dims_token + pad_dims_template,
            ),
            template_distances=F.pad(
                self.template_distances,
                pad=pad_dims_token + pad_dims_token + pad_dims_template,
            ),
            template_unit_vector=F.pad(
                self.template_unit_vector,
                # This field has a final dimension of size 3, which we shouldn't pad
                pad=(0, 0) + pad_dims_token + pad_dims_token + pad_dims_template,
            ),
        )

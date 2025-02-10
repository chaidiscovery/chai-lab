# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator

import torch
from torch import Tensor
from torch.nn import functional as F

from chai_lab.data import residue_constants as rc
from chai_lab.data.dataset.structure.all_atom_residue_tokenizer import (
    AllAtomResidueTokenizer,
)
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.dataset.templates.align import align_1d, align_2d
from chai_lab.data.dataset.templates.load import LoadedTemplate, get_template_data
from chai_lab.data.parsing.msas.aligned_pqt import hash_sequence
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.parsing.templates.m8 import parse_m8_to_template_hits
from chai_lab.data.parsing.templates.template_hit import TemplateHit
from chai_lab.data.sources.rdkit import RefConformerGenerator
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

    @classmethod
    def merge(
        cls,
        templates: list["TemplateContext"],
    ) -> "TemplateContext":
        """Merge template contexts along the template dimensions."""
        # n_token can be simply concatenated
        logger.debug(f"Merging {len(templates)} templates")

        # Handle case where we get an empty list (no templates to merge)
        if len(templates) == 0:
            return cls.empty(n_templates=4, n_tokens=1)

        # Pad each template_restype's template_dimension to match the largest
        # NOTE count num_templates here, NOT num_nonnull_templates
        n_templates_new: int = max(t.num_templates for t in templates)
        padded_templates = [t.pad(max_templates=n_templates_new) for t in templates]
        new_template_restype = torch.cat(
            [t.template_restype for t in padded_templates],
            dim=1,  # Concat on sequence dim
        )
        new_template_pseudo_beta_mask = torch.cat(
            [t.template_pseudo_beta_mask for t in padded_templates],
            dim=1,
        )
        new_template_backbone_frame_mask = torch.cat(
            [t.template_backbone_frame_mask for t in padded_templates],
            dim=1,
        )

        # Number of tokens after concatenation along token dim
        n_token_new = new_template_restype.shape[1]

        # n_token x n_token must be tiled into a square matrix
        # These indices like [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 3, ...] indicate the region
        # of the square matrix that corresponds to each template.
        template_indices = torch.repeat_interleave(
            input=torch.arange(len(templates), device=new_template_restype.device),
            repeats=torch.tensor([t.template_restype.shape[-1] for t in templates]),
        )
        # Sample template and token dim
        assert template_indices.shape[0] == n_token_new

        new_template_distances = torch.zeros(
            n_templates_new, n_token_new, n_token_new, dtype=torch.float32
        )
        new_template_unit_vector = torch.zeros(
            n_templates_new, n_token_new, n_token_new, 3, dtype=torch.float32
        )

        # For each template, find the block that it corresponds to and copy in the data
        for i, t in enumerate(templates):
            m = template_indices == i
            mask = m[:, None] * m[None, :]
            idx = torch.arange(t.template_distances.shape[0])
            new_template_distances[idx.unsqueeze(1), mask] = (
                t.template_distances.flatten(1, 2)
            )
            new_template_unit_vector[idx.unsqueeze(1), mask] = (
                t.template_unit_vector.flatten(1, 2)
            )

        return cls(
            template_restype=new_template_restype,
            template_pseudo_beta_mask=new_template_pseudo_beta_mask,
            template_backbone_frame_mask=new_template_backbone_frame_mask,
            template_distances=new_template_distances,
            template_unit_vector=new_template_unit_vector,
        )

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

    @classmethod
    def from_loaded_templates(
        cls,
        n_tokens: int,
        templates: list[LoadedTemplate],
        pad_to_n_templates: int = 0,
        apply_crop: bool = False,
    ) -> "TemplateContext":
        """Aligns templates into a single tensor whose size matches the query.

        Setting apply_crop controls whether the output is aligned wrt the cropped query
        or the full query.

        If pad_to_n_templates is set to 0, then we do not pad the templates, and just
        return whatever number of templates are input.
        """
        if len(templates) == 0:
            # If pad to n templates is 0, then we return one single dummy template
            # If pad_to_n_templates is nonzero, then we return that many dummy templates
            return cls.empty(n_templates=max(1, pad_to_n_templates), n_tokens=n_tokens)

        if pad_to_n_templates == 0:
            # If padding is disabled, simply set to the number of input templates
            pad_to_n_templates = len(templates)
        else:
            # If explicitly padding to a set number of templates, check that we haven't
            # provided too many templates.
            assert len(templates) <= pad_to_n_templates
        assert pad_to_n_templates > 0

        # Select the indices
        alignment_indices = [
            t.cropped_template_query_match_indices
            if apply_crop
            else t.template_query_match_indices
            for t in templates
        ]
        if (
            max_align_idx := max([torch.max(a).item() for a in alignment_indices])
        ) >= n_tokens:
            raise ValueError(
                "Template alignment indices exceed number of tokens: "
                f"{apply_crop=} {max_align_idx=} >= {n_tokens=}."
            )

        template_restype = align_1d(
            n_tokens,
            tensors=[t.template_restype for t in templates],
            indices=alignment_indices,
            pad_to_n_templates=pad_to_n_templates,
            default_dtype=torch.int32,
            fill_value=rc.residue_types_with_nucleotides_order["-"],
        )
        template_pseudo_beta_mask = align_1d(
            n_tokens,
            tensors=[t.template_pseudo_beta_mask for t in templates],
            indices=alignment_indices,
            pad_to_n_templates=pad_to_n_templates,
            default_dtype=torch.bool,
        )
        template_backbone_frame_mask = align_1d(
            n_tokens,
            tensors=[t.template_backbone_frame_mask for t in templates],
            indices=alignment_indices,
            pad_to_n_templates=pad_to_n_templates,
            default_dtype=torch.bool,
        )
        template_distances = align_2d(
            n_tokens,
            tensors=[t.template_pseudo_beta_distances for t in templates],
            indices=alignment_indices,
            pad_to_n_templates=pad_to_n_templates,
            default_dtype=torch.float,
        )
        template_unit_vector = align_2d(
            n_tokens,
            tensors=[t.template_unit_vector for t in templates],
            indices=alignment_indices,
            pad_to_n_templates=pad_to_n_templates,
            default_dtype=torch.float,
            addtl_dims=[3],  # Unit vectors are 3 values per entry
        )
        return TemplateContext(
            template_restype=template_restype,
            template_pseudo_beta_mask=template_pseudo_beta_mask,
            template_backbone_frame_mask=template_backbone_frame_mask,
            template_distances=template_distances,
            template_unit_vector=template_unit_vector,
        )


def get_template_context(
    chains: list[Chain],
    template_hits_m8: Path,
    use_sequence_hash_for_lookup: bool = False,
    template_cif_cache_folder: Path | None = None,
) -> TemplateContext:
    """
    For each example, loads templates for cropped chain, collate the templates.

    Uses the following logic:
    - Empty crops get dropped
    - Crops that don't have templates (no hits, non-protein entity, etc.) have a empty
      template loaded to preserve alignment when merging templates.
    - Crops with templates are loaded from a3m file, aligned back to the query, cropped
      down, and blown out based on which tokens may have been tokenized per atom.
    """
    templates: list[TemplateContext] = []
    tok = AllAtomResidueTokenizer(RefConformerGenerator())

    for chain in chains:
        # Below, any time we append a template context, we must ensure that the
        # size matches the cropped original context such that the original
        # contexts and the templates line up.
        structure_context: AllAtomStructureContext = chain.structure_context

        _, token_residue_index_zeroed = torch.unique(
            structure_context.token_residue_index,
            return_inverse=True,
        )

        # Get the template hits for this chain
        # Create a new iterator of template hits using the arguments such that we never
        # consume an already-advanced iterator. This has a slight cost of loading in the
        # a3m file every time, but this is relatively small compared to the cost of
        # actually aligning the templates to the query, which we can do lazily with an
        # iterator.
        loaded_templates: list[LoadedTemplate] = []
        if chain.entity_data.entity_type == EntityType.PROTEIN:
            # Create an iterator over tempalte hits
            template_hits: Iterator[TemplateHit] = parse_m8_to_template_hits(
                (
                    chain.entity_data.entity_name
                    if not use_sequence_hash_for_lookup
                    else hash_sequence(chain.entity_data.sequence)
                ),
                chain.entity_data.sequence,
                template_hits_m8,
                template_cif_folder=template_cif_cache_folder,
            )
            # Load the template data
            loaded_templates = get_template_data(
                template_hits=template_hits,
                query_crop_indices=torch.arange(structure_context.num_tokens),
                tokenizer=tok,
                strict_subsequence_check=True,
                drop_unresolved_from_hits=True,
            )

        # This can get triggered if it is a non-protein or if it's a protein with no
        # template hits that pass filters.
        if len(loaded_templates) == 0:
            # No templates - add an empty TemplateContext of matching size
            empty = TemplateContext.empty(
                n_templates=1, n_tokens=structure_context.num_tokens
            ).index_select(token_residue_index_zeroed)
            templates.append(empty)
            continue

        # Align templates into a TemplateContext whose size matches UNCROPPED query.
        template_context = TemplateContext.from_loaded_templates(
            n_tokens=structure_context.num_tokens,
            pad_to_n_templates=0,
            templates=loaded_templates,
            apply_crop=False,
        )

        # Blow-out to match cases of per-atom tokenization
        blown_out_template_context = template_context.index_select(
            token_residue_index_zeroed
        )
        templates.append(blown_out_template_context)

    for chain, template_ctx in zip(chains, templates, strict=True):
        assert chain.num_tokens == template_ctx.num_tokens
    # Merge templates across the chains
    return TemplateContext.merge(templates)

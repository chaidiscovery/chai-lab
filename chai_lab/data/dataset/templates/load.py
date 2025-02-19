# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
"""
Consider the following example:

Crop indices:               |-------------------------------------|
Query: ------------------------------------------------------------------------
Hit on 1B2D:                       |--------------------------|
Hit start                          |
Hit end                                                       |

In the LoadedTemplate structure:
- template_hit_indices             |                          |  (wrt full hit sequence)
- cropped_query_match_indices      |                          |  (wrt cropped query)
- cropped_query_match_mask  -------*************__*************---- (wrt cropped query)
"""

import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterator

import gemmi
import torch
from aiohttp import ClientResponseError
from torch import Tensor
from urllib3.exceptions import MaxRetryError, NameResolutionError

from chai_lab.data import residue_constants as rc
from chai_lab.data.dataset.structure.all_atom_residue_tokenizer import (
    AllAtomResidueTokenizer,
)
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.features.token_utils import (
    get_centre_positions_and_mask,
    get_token_reference_atom_positions_and_mask,
)
from chai_lab.data.io.rcsb import download_cif_file
from chai_lab.data.parsing.structure.all_atom_entity_data import (
    AllAtomEntityData,
    structure_to_entities_data,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.parsing.templates.template_hit import TemplateHit
from chai_lab.tools.rigid import Rigid
from chai_lab.utils.tensor_utils import cdist
from chai_lab.utils.typing import Bool, Float, Int, typecheck

logger = logging.getLogger(__name__)
logger.setLevel(level=logging.INFO)


@typecheck
@dataclass(frozen=True)
class LoadedTemplate:
    """Template data structure.

    Compared to the TemplateHit structure under hhr, which holds hits discovered by
    HHSearch, this structure holds the actual templates that we want to featurize. To
    featurize a template, we are interested in its:
    * sequence
    * N-Ca-C coordinates (backbone) + Cb coordinate
    * Bacbone mask indicating existence of Cb
    * Backbone frame mask indicating whether all atom exist for backbone frame
    """

    query_crop_indices: Int[Tensor, "c"]  # != region of template match
    template_hit: TemplateHit
    template_hit_structure_context: AllAtomStructureContext  # Not cropped
    # Cropping/indexing into the above is done by template_hit.indices_hit

    def __post_init__(self):
        # Query crop indices should be strictly increasing
        diffs = self.query_crop_indices[1:] - self.query_crop_indices[:-1]
        assert torch.all(
            diffs > 0
        ), f"Query crop indices not strictly increasing: {self.query_crop_indices}"

        # Template hit is fully contained within query indices
        query_range = set(self.query_crop_indices.tolist())
        hit_range = set(self.template_hit.indices_query.tolist())
        assert hit_range.intersection(
            query_range
        ), f"Hit {self.template_hit.query_start_end} does not overlap {self.query_identifier} crop indices"

    @property
    def query_identifier(self) -> str:
        return self.template_hit.query_pdb_id

    @property
    def hit_identifier(self) -> str:
        return f"{self.template_hit.pdb_id}|{self.template_hit.chain_id}"

    @property
    def template_hit_indices(self) -> Int[Tensor, "n_tokens"]:
        """Indices of hit within full hit sequence."""
        return self.template_hit.indices_hit

    @property
    def template_query_match_indices(self) -> Int[Tensor, "n_tokens"]:
        """Indices of query residues corresponding to hit."""
        return self.template_hit.indices_query

    @property
    def cropped_template_query_match_indices(self) -> Int[Tensor, "n_tokens"]:
        """Indices of query residues corresponding to hit within cropped template."""
        retval = torch.searchsorted(
            self.query_crop_indices, self.template_hit.indices_query
        )
        assert torch.all(retval >= 0)
        return retval

    # The below are the "features" that are eventually given to the model as input

    @property
    def template_restype(self) -> Int[Tensor, "n_tokens"]:
        # Get the residue type and crop it
        raw_restype = self.template_hit_structure_context.token_residue_type.squeeze(0)
        retval = raw_restype[self.template_hit_indices]
        retval[~self.template_hit.hit_valid_mask] = (
            rc.residue_types_with_nucleotides_order["-"]
        )
        return retval

    @property
    def template_pseudo_beta_mask(self) -> Bool[Tensor, "n_tokens"]:
        _, mask = get_token_reference_atom_positions_and_mask(
            atom_pos=self.template_hit_structure_context.atom_gt_coords.unsqueeze(0),
            atom_mask=self.template_hit_structure_context.atom_exists_mask.unsqueeze(0),
            token_reference_atom_index=self.template_hit_structure_context.token_ref_atom_index.unsqueeze(
                0
            ),
            token_exists_mask=self.template_hit_structure_context.token_exists_mask.unsqueeze(
                0
            ),
        )
        retval = mask.squeeze(0)[self.template_hit_indices]
        # hit_valid_mask is True if valid; un-set things that are not valid
        return retval & self.template_hit.hit_valid_mask

    @property
    def template_pseudo_beta_distances(self) -> Bool[Tensor, "n_tokens n_tokens"]:
        coords, mask = get_token_reference_atom_positions_and_mask(
            atom_pos=self.template_hit_structure_context.atom_gt_coords.unsqueeze(0),
            atom_mask=self.template_hit_structure_context.atom_exists_mask.unsqueeze(0),
            token_reference_atom_index=self.template_hit_structure_context.token_ref_atom_index.unsqueeze(
                0
            ),
            token_exists_mask=self.template_hit_structure_context.token_exists_mask.unsqueeze(
                0
            ),
        )

        coords = coords.squeeze(0)[self.template_hit_indices]
        mask = (
            mask.squeeze(0)[self.template_hit_indices]  # Reference positions mask
            & self.template_hit.hit_valid_mask  # Hit valid mask (not a gap)
        )
        # Replace masked regions with some large value such that when we one-hot encode
        # the distance matrix, the masked values land in the last "everything too far"
        # block. The default one-hot encoding caps at 50.75, so the masked value of 100
        # is well above this.
        distances = cdist(coords.unsqueeze(0)).squeeze(0)
        mask_block = mask[:, None] * mask[None, :]
        distances.masked_fill_(~mask_block, 100.0)
        return distances

    @property
    def template_backbone_frame_mask(self) -> Bool[Tensor, "n_tokens"]:
        retval = self.template_hit_structure_context.token_backbone_frame_mask.squeeze(
            0
        )
        return retval[self.template_hit_indices] & self.template_hit.hit_valid_mask

    @property
    def template_unit_vector(self) -> Float[Tensor, "n_tokens n_tokens 3"]:
        eps = 1e-12
        template_mask_2d = (
            self.template_backbone_frame_mask[..., None]
            * self.template_backbone_frame_mask[..., None, :]
        )

        # Start by compuing the coordinates for all backbone frames
        backbone_frame_indices = (
            self.template_hit_structure_context.token_backbone_frame_index
        )
        n, ca, c = [
            i.squeeze(-1) for i in torch.split(backbone_frame_indices, 1, dim=-1)
        ]
        # These are per-token wrt the original full context
        backbone_positions_and_mask = [
            get_token_reference_atom_positions_and_mask(
                atom_pos=self.template_hit_structure_context.atom_gt_coords.unsqueeze(
                    0
                ),
                atom_mask=self.template_hit_structure_context.atom_exists_mask.unsqueeze(
                    0
                ),
                token_reference_atom_index=idx.unsqueeze(0),
                token_exists_mask=self.template_hit_structure_context.token_backbone_frame_mask.unsqueeze(
                    0
                ),
            )
            for idx in [n, ca, c]
        ]
        # Crop both the positions and the masks to the regions of the template hit
        # NOTE we don't use the masks generated from the above because we mask separately below
        n_pos, ca_pos, c_pos = [
            pos[:, self.template_hit_indices]
            for pos, _mask in backbone_positions_and_mask
        ]

        rigids = Rigid.make_transform_from_reference(
            n_xyz=n_pos,
            ca_xyz=ca_pos,
            c_xyz=c_pos,
            eps=eps,
        )
        points = rigids.get_trans()
        rigid_vec = rigids[..., None].invert_apply(points)

        inv_distance_scalar = torch.rsqrt(eps + torch.sum(rigid_vec**2, dim=-1))

        inv_distance_scalar = inv_distance_scalar * template_mask_2d
        unit_vector = rigid_vec * inv_distance_scalar[..., None]
        return unit_vector.squeeze(0)


def _get_entity_data(
    pdb_id_or_path: str | Path,
    chain_id: str,
    subchain_id: str | None = None,
) -> AllAtomEntityData:
    if not Path(pdb_id_or_path).is_file():
        with TemporaryDirectory() as tmpdir:
            local_path = download_cif_file(str(pdb_id_or_path), Path(tmpdir))
            gemmi_structure = gemmi.read_structure(str(local_path))
    else:
        gemmi_structure = gemmi.read_structure(str(pdb_id_or_path))

    entities_data = structure_to_entities_data(
        structure=gemmi_structure,
        chain_ids=[chain_id],
        subchain_ids=[subchain_id] if subchain_id is not None else None,
        make_assembly=False,
    )
    # Protein only
    entities_data = [e for e in entities_data if e.entity_type == EntityType.PROTEIN]
    if not len(entities_data) == 1:
        raise ValueError(
            f"Expected exactly one protein entity for {pdb_id_or_path} chain {chain_id}, but got {len(entities_data)}"
        )
    return entities_data[0]


@typecheck
def get_template_data(
    template_hits: Iterator[TemplateHit],
    tokenizer: AllAtomResidueTokenizer,
    query_crop_indices: Int[Tensor, "tok"],
    max_loaded_templates: int = 4,
    strict_subsequence_check: bool = True,
    fully_contained_only: bool = False,
    drop_unresolved_from_hits: bool = True,
) -> list[LoadedTemplate]:
    """Get template data for the given query PDB chain subject to the date cutoff.

    Setting max_loaded_templates to be equal to subsample_templates will result in
    directly loading max_loaded_templates and returning them without subsampling.

    subsample_randomly controls whether we subsample from max_loaded_templates to
    subsample_templates randomly (True) or by taking the first few templates (False)

    If strict_subsequence_check is True, we check that the template sequence is a
    subsequence of the template hit sequence and skip instances where it is not.

    If fully_contained_only is True, then we only keep template hits that are fully
    contained within the crop indices; otherwise we require nonzero overlap. Enabling
    this applies a stricter thresholding to templates and removes any templates that
    have overhangs outside of the cropped region.
    """
    query_crop_indices_set = set(query_crop_indices.tolist())

    # Gather template hits
    template_data: list[LoadedTemplate] = []
    dropped_hits: dict[str, list[TemplateHit]] = defaultdict(list)

    # If the iterator is empty, we won't run this loop
    for template_hit in template_hits:
        # Keep only hits fully contained within query indices
        hit_query_indices = {i for i in template_hit.indices_query.tolist() if i >= 0}
        if fully_contained_only:
            overlaps = hit_query_indices.issubset(query_crop_indices_set)
        else:
            overlaps = len(hit_query_indices.intersection(query_crop_indices_set)) > 0
        if not overlaps:
            logger.debug(f"{template_hit} does not intersect crop indices")
            continue

        # Load the entity data for the hit
        try:
            template_entity_data = _get_entity_data(
                pdb_id_or_path=(
                    template_hit.pdb_id
                    if template_hit.cif_path is None
                    else template_hit.cif_path
                ),
                chain_id=template_hit.chain_id,  # Templates are referenced by chain, NOT subchain
            )
        except (ValueError, ClientResponseError, NameResolutionError, MaxRetryError):
            logger.info(f"Failed to load entity data for {template_hit}")
            continue
        except Exception:
            logger.exception(
                f"[Unknown case - you should debug this! - but we'll handle gracefully] Failed to load entity data for {template_hit}"
            )
            continue

        if template_entity_data.has_modifications:
            # NOTE we drop hits with modifications because these get tokenized per-atom,
            # so in order to get them to match up with the original query/input
            # sequence, we would need a way to collapse the per-token tokenization into
            # a single token.
            logger.debug(f"Skipping {template_hit} as it has modified residues")
            continue

        template_structure_context = tokenizer.tokenize_entity(template_entity_data)
        if template_structure_context is None:
            logger.warning(f"Skipping {template_hit} as it failed to tokenize")
            continue

        # Template hits are
        if drop_unresolved_from_hits:
            _, mask = get_centre_positions_and_mask(
                atom_gt_coords=template_structure_context.atom_gt_coords,
                atom_exists_mask=template_structure_context.atom_exists_mask,
                token_centre_atom_index=template_structure_context.token_centre_atom_index,
                token_exists_mask=template_structure_context.token_exists_mask,
            )
            (resolved_indices,) = torch.where(mask)
            template_structure_context = template_structure_context.index_select(
                resolved_indices
            )

        # N.B. it is possible that entries tokenized by atoms cause issues here due to
        # the crop being small ut the number of tokens that come through being large.
        # ALternatively, this might be because the chain is bigger than the subchain
        # so... maybe we have to work on subchain level?
        template = LoadedTemplate(
            query_crop_indices=query_crop_indices,
            template_hit=template_hit,
            template_hit_structure_context=template_structure_context,
        )

        # Check that the loaded version based on an AllAtomStructureContext matches what
        # we expect from the "raw" TemplateHit
        if strict_subsequence_check:
            try:
                original_alignment = "".join(
                    [
                        c
                        for c in template_hit.query_seq_realigned
                        if c.isupper() or not c.isalpha()  # Keep upper and "-"
                    ]
                )
                context_alignment = "".join(
                    [
                        rc.residue_types_with_nucleotides[i]
                        for i in template.template_restype.tolist()
                    ]
                )
                # NOTE before these are all loaded into a TemplateContext, they do not
                # contain any flanking gap chars "-" so we do not check that the flanking
                # gaps are equal.
                if context_alignment not in original_alignment:
                    logger.warning(
                        f"Skipping {template_hit} due to mismatched sequences: {context_alignment=} {original_alignment=}"
                    )
                    continue
            except IndexError:
                logger.warning(
                    f"Skipping {template_hit} due to exception when checking sequences",
                    exc_info=True,
                )
                continue

        # Add the hit to the list of hits, stopping once we have enough
        template_data.append(template)
        if len(template_data) >= max_loaded_templates:
            break

    if len(dropped_hits) > 0:
        # If we have any dropped hits, that means we have encountered at least one
        # template hit so this will be defined.
        assert isinstance(template_hit, TemplateHit)
        pdb_id = template_hit.query_pdb_id
        drop_count = {k: len(v) for k, v in dropped_hits.items()}
        logger.debug(
            f"Templates for {pdb_id} | {len(template_data)} remain, dropped hits: {drop_count}"
        )

    logger.info(
        f"Loaded {len(template_data)} templates: {[t.hit_identifier for t in template_data]}"
    )
    return template_data

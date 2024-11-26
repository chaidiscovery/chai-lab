# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import dataclass
from typing import Any, Final

from chai_lab.data.dataset.constraints.restraint_context import RestraintContext
from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.dataset.templates.context import TemplateContext

logger = logging.getLogger(__name__)

MAX_MSA_DEPTH: Final[int] = 16_384
MAX_NUM_TEMPLATES: Final[int] = 4


@dataclass
class AllAtomFeatureContext:
    """
    Feature contexts are produced by datasets. Multiple feature contexts are passed to
    collator, which transforms them into a batch (by padding and stacking them).
    """

    # Metadata: these are not padded and batched
    chains: list[Chain]
    # Contexts: these are what get padded and batched
    structure_context: AllAtomStructureContext
    msa_context: MSAContext
    profile_msa_context: MSAContext
    template_context: TemplateContext
    embedding_context: EmbeddingContext | None
    restraint_context: RestraintContext

    def __str__(self) -> str:
        chains_info = [str(chain) for chain in self.chains]
        return f"{self.__class__.__name__}(chains={chains_info})"

    def pad(
        self,
        n_tokens: int,
        n_atoms: int,
    ) -> "AllAtomFeatureContext":
        return AllAtomFeatureContext(
            # Metadata
            chains=self.chains,
            # Contexts
            structure_context=self.structure_context.pad(
                n_tokens=n_tokens,
                n_atoms=n_atoms,
            ),
            msa_context=self.msa_context.pad(
                max_num_tokens=n_tokens,
                max_msa_depth=MAX_MSA_DEPTH,
            ),
            profile_msa_context=self.profile_msa_context.pad(
                max_num_tokens=n_tokens,
                # max_msa_depth=MAX_MSA_DEPTH,
            ),
            template_context=self.template_context.pad(
                max_tokens=n_tokens,
                max_templates=MAX_NUM_TEMPLATES,
            ),
            embedding_context=(
                self.embedding_context.pad(max_tokens=n_tokens)
                if self.embedding_context is not None
                else None
            ),
            restraint_context=self.restraint_context.pad(max_tokens=n_tokens),
        )

    def to_dict(self) -> dict[str, Any]:
        msa_context_dict = dict(
            msa_tokens=self.msa_context.tokens,
            msa_mask=self.msa_context.mask,
            msa_deletion_matrix=self.msa_context.deletion_matrix,
            msa_pairkey=self.msa_context.pairing_key_hash,
            msa_sequence_source=self.msa_context.sequence_source,
            main_msa_tokens=self.profile_msa_context.tokens,
            main_msa_mask=self.profile_msa_context.mask,
            main_msa_deletion_matrix=self.profile_msa_context.deletion_matrix,
        )
        return {
            **self.structure_context.to_dict(),
            **msa_context_dict,
            **self.template_context.to_dict(),
            **(self.embedding_context.to_dict() if self.embedding_context else {}),
            **self.restraint_context.to_dict(),
        }

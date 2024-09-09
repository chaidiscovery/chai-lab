from dataclasses import dataclass

import torch
from einops import rearrange, repeat
from torch import Tensor

from chai_lab.data.parsing.msas.data_source import (
    MSADataSource,
    msa_dataset_source_to_int,
)
from chai_lab.data.parsing.msas.species import UNKNOWN_SPECIES
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order
from chai_lab.utils.defaults import default
from chai_lab.utils.typing import Bool, Int32, UInt8, typecheck


@typecheck
@dataclass
class MSAContext:
    # MSA-level
    dataset_source: MSADataSource

    # token level
    tokens: UInt8[Tensor, "msa_depth n_tokens"]
    species: Int32[Tensor, "msa_depth n_tokens"]
    deletion_matrix: UInt8[Tensor, "msa_depth n_tokens"]
    mask: Bool[Tensor, "msa_depth n_tokens"]
    sequence_source: UInt8[Tensor, "msa_depth n_tokens"]
    is_paired_mask: Bool[Tensor, "msa_depth"]

    @property
    def depth(self) -> int:
        depth, _ = self._dims
        return depth

    @property
    def num_tokens(self) -> int:
        _, num_tokens = self._dims
        return num_tokens

    @property
    def _dims(self) -> torch.Size:
        return self.tokens.shape

    @property
    def paired_msa_depth(self) -> Int32[Tensor, "b"]:
        return (self.mask.any(dim=-1) & self.is_paired_mask).sum(dim=-1)

    def __getitem__(self, subscript: tuple) -> "MSAContext":
        # enforce typing on item
        if not (
            isinstance(subscript, tuple)
            and ((len(subscript) == 2) or subscript[0] is Ellipsis)
        ):
            raise TypeError(
                "Subscript must be a tuple with 2 elements or have an ellipsis."
            )

        is_paired_mask = repeat(
            self.is_paired_mask,
            "msa_depth -> msa_depth n_tokens",
            n_tokens=self.num_tokens,
        )
        return MSAContext(
            dataset_source=self.dataset_source,
            tokens=self.tokens[subscript],
            species=self.species[subscript],
            deletion_matrix=self.deletion_matrix[subscript],
            sequence_source=self.sequence_source[subscript],
            mask=self.mask[subscript],
            is_paired_mask=is_paired_mask[subscript].any(dim=-1),
        )

    def pad(
        self,
        max_num_tokens: int | None = None,
        max_msa_depth: int | None = None,
    ) -> "MSAContext":
        max_num_tokens = default(max_num_tokens, self.num_tokens)
        assert self.num_tokens <= max_num_tokens

        max_msa_depth = default(max_msa_depth, self.depth)
        assert self.depth <= max_msa_depth

        pad_dims = (0, max_num_tokens - self.num_tokens, 0, max_msa_depth - self.depth)
        return MSAContext(
            dataset_source=self.dataset_source,
            tokens=torch.nn.functional.pad(
                self.tokens,
                pad_dims,
                value=residue_types_with_nucleotides_order[":"],
            ),
            species=torch.nn.functional.pad(
                self.species,
                pad_dims,
                value=UNKNOWN_SPECIES,
            ),
            deletion_matrix=torch.nn.functional.pad(
                self.deletion_matrix,
                pad_dims,
                value=0,  # No deletions
            ),
            mask=torch.nn.functional.pad(
                self.mask,
                pad_dims,
                value=False,
            ),
            sequence_source=torch.nn.functional.pad(
                self.sequence_source,
                pad_dims,
                value=msa_dataset_source_to_int[MSADataSource.NONE],
            ),
            is_paired_mask=torch.nn.functional.pad(
                self.is_paired_mask,
                (0, max_msa_depth - self.depth),
                value=False,
            ),
        )

    @typecheck
    def apply_mask(self, mask: Bool[Tensor, "msa_depth n_tokens"]) -> "MSAContext":
        return MSAContext(
            dataset_source=self.dataset_source,
            tokens=self.tokens.masked_fill(
                ~mask, residue_types_with_nucleotides_order[":"]
            ),
            species=self.species.masked_fill(~mask, UNKNOWN_SPECIES),
            deletion_matrix=self.deletion_matrix.masked_fill(~mask, 0),
            mask=self.mask.masked_fill(~mask, False),
            sequence_source=self.sequence_source.masked_fill(
                ~mask, value=msa_dataset_source_to_int[MSADataSource.NONE]
            ),
            is_paired_mask=self.is_paired_mask.masked_fill(~mask.any(dim=-1), False),
        )

    @classmethod
    def cat(
        cls,
        msas: list["MSAContext"],
        dataset_source: MSADataSource | None = None,
        dim=-1,
    ) -> "MSAContext":
        if dataset_source is None:
            dataset_sources = set([msa.dataset_source for msa in msas])
            assert len(dataset_sources) == 1 or dataset_sources == {
                MSADataSource.MAIN,
                MSADataSource.PAIRED,
            }, "all MSAs must have the same datasource or be MAIN and PAIRED"
            dataset_source = dataset_sources.pop()

        assert dim == -1 or dim >= 0, "dim < 0 not implemented except for -1"
        if 0 <= dim < 1:
            is_paired_mask = torch.cat([msa.is_paired_mask for msa in msas], dim=dim)
        else:
            assert len(msas) > 0
            is_paired_mask = msas[0].is_paired_mask

        return MSAContext(
            dataset_source=dataset_source,
            tokens=torch.cat([msa.tokens for msa in msas], dim=dim),
            species=torch.cat([msa.species for msa in msas], dim=dim),
            deletion_matrix=torch.cat([msa.deletion_matrix for msa in msas], dim=dim),
            sequence_source=torch.cat([msa.sequence_source for msa in msas], dim=dim),
            mask=torch.cat([msa.mask for msa in msas], dim=dim),
            is_paired_mask=is_paired_mask,
        )

    @classmethod
    @typecheck
    def create(
        cls,
        dataset_source: MSADataSource,
        tokens: UInt8[Tensor, "n_tokens"],
    ) -> "MSAContext":
        """
        Creates an MSA comprised of a single sequence.
        """
        tokens_for_msa = rearrange(tokens, "n_tokens -> 1 n_tokens")
        return MSAContext(
            dataset_source=dataset_source,
            tokens=tokens_for_msa,
            species=torch.full_like(tokens_for_msa, UNKNOWN_SPECIES, dtype=torch.int32),
            deletion_matrix=torch.zeros_like(tokens_for_msa, dtype=torch.uint8),
            mask=torch.ones_like(tokens_for_msa, dtype=torch.bool),
            sequence_source=torch.full_like(
                tokens_for_msa,
                fill_value=msa_dataset_source_to_int[dataset_source],
            ),
            is_paired_mask=torch.zeros((1,), dtype=torch.bool),
        )

    @classmethod
    def create_empty(cls, n_tokens: int, depth: int = 0) -> "MSAContext":
        dims = (depth, n_tokens)
        return MSAContext(
            dataset_source=MSADataSource.NONE,
            tokens=torch.full(
                dims, residue_types_with_nucleotides_order[":"], dtype=torch.uint8
            ),
            species=torch.full(dims, UNKNOWN_SPECIES, dtype=torch.int32),
            deletion_matrix=torch.zeros(dims, dtype=torch.uint8),  # No deletions
            mask=torch.zeros(dims, dtype=torch.bool),
            sequence_source=torch.full(
                dims,
                fill_value=msa_dataset_source_to_int[MSADataSource.NONE],
                dtype=torch.uint8,
            ),
            is_paired_mask=torch.zeros((depth,), dtype=torch.bool),
        )

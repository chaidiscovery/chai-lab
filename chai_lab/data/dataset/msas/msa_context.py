# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.parsing.msas.data_source import (
    MSADataSource,
    msa_dataset_source_to_int,
)
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order
from chai_lab.utils.defaults import default
from chai_lab.utils.typing import Bool, Int32, UInt8, typecheck

NO_PAIRING_KEY = -999991


@typecheck
@dataclass
class MSAContext:
    # token level
    tokens: UInt8[Tensor, "msa_depth n_tokens"]
    pairing_key_hash: Int32[Tensor, "msa_depth n_tokens"]
    deletion_matrix: UInt8[Tensor, "msa_depth n_tokens"]
    mask: Bool[Tensor, "msa_depth n_tokens"]
    sequence_source: UInt8[Tensor, "msa_depth n_tokens"]

    @property
    def depth(self) -> int:
        depth, _ = self.tokens.shape
        return depth

    @property
    def num_tokens(self) -> int:
        _, num_tokens = self.tokens.shape
        return num_tokens

    def __getitem__(self, subscript: tuple) -> "MSAContext":
        assert (
            isinstance(subscript, tuple) and len(subscript) == 2
        ), "Subscript must be a tuple with 2 elements or have an ellipsis."

        return MSAContext(
            tokens=self.tokens[subscript],
            pairing_key_hash=self.pairing_key_hash[subscript],
            deletion_matrix=self.deletion_matrix[subscript],
            sequence_source=self.sequence_source[subscript],
            mask=self.mask[subscript],
        )

    def take_rows_with_padding(self, row_indices_with_nones: list[int | None]):
        """
        allows specifying index=None, which will be filled with empty sequence,
        helpful to align multiple sequences
        """
        padded = self.pad(max_msa_depth=self.depth + 1)  # add one empty row, idx=-1
        filled_indices = torch.asarray(
            [-1 if x is None else x for x in row_indices_with_nones], dtype=torch.long
        )  # idx=-1 for new sequence
        return padded[filled_indices, :]

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

        def pad2d(tensor: torch.Tensor, pad_value) -> torch.Tensor:
            return torch.nn.functional.pad(tensor, pad_dims, value=pad_value)

        return MSAContext(
            tokens=pad2d(self.tokens, residue_types_with_nucleotides_order[":"]),
            pairing_key_hash=pad2d(self.pairing_key_hash, pad_value=NO_PAIRING_KEY),
            deletion_matrix=pad2d(self.deletion_matrix, pad_value=0),  # No deletions
            mask=pad2d(self.mask, pad_value=False),
            sequence_source=pad2d(
                self.sequence_source,
                msa_dataset_source_to_int[MSADataSource.NONE],
            ),
        )

    @typecheck
    @classmethod
    def cat(cls, msas: list["MSAContext"], dim: int) -> "MSAContext":
        assert dim in (0, 1, -1), f"unexpected {dim=}, not accepting dim < -1"

        return MSAContext(
            tokens=torch.cat([msa.tokens for msa in msas], dim=dim),
            pairing_key_hash=torch.cat([msa.pairing_key_hash for msa in msas], dim=dim),
            deletion_matrix=torch.cat([msa.deletion_matrix for msa in msas], dim=dim),
            sequence_source=torch.cat([msa.sequence_source for msa in msas], dim=dim),
            mask=torch.cat([msa.mask for msa in msas], dim=dim),
        )

    @typecheck
    def apply_mask(self, mask: Bool[Tensor, "msa_depth n_tokens"]) -> "MSAContext":
        return MSAContext(
            tokens=self.tokens.masked_fill(
                ~mask, residue_types_with_nucleotides_order[":"]
            ),
            pairing_key_hash=self.pairing_key_hash.masked_fill(~mask, NO_PAIRING_KEY),
            deletion_matrix=self.deletion_matrix.masked_fill(~mask, 0),
            mask=self.mask.masked_fill(~mask, False),
            sequence_source=self.sequence_source.masked_fill(
                ~mask, value=msa_dataset_source_to_int[MSADataSource.NONE]
            ),
        )

    @classmethod
    @typecheck
    def create_single_seq(
        cls,
        dataset_source: MSADataSource,
        tokens: UInt8[Tensor, "n_tokens"],
    ) -> "MSAContext":
        """
        Creates an MSA comprised of a single sequence.
        """
        tokens_for_msa = rearrange(tokens, "n_tokens -> 1 n_tokens")
        return MSAContext(
            tokens=tokens_for_msa,
            pairing_key_hash=torch.full_like(
                tokens_for_msa, NO_PAIRING_KEY, dtype=torch.int32
            ),
            deletion_matrix=torch.zeros_like(tokens_for_msa, dtype=torch.uint8),
            mask=torch.ones_like(tokens_for_msa, dtype=torch.bool),
            sequence_source=torch.full_like(
                tokens_for_msa,
                fill_value=msa_dataset_source_to_int[dataset_source],
            ),
        )

    @classmethod
    def create_empty(cls, n_tokens: int, depth: int = 0) -> "MSAContext":
        dims = (depth, n_tokens)
        return MSAContext(
            tokens=torch.full(
                dims, residue_types_with_nucleotides_order[":"], dtype=torch.uint8
            ),
            pairing_key_hash=torch.full(dims, NO_PAIRING_KEY, dtype=torch.int32),
            deletion_matrix=torch.zeros(dims, dtype=torch.uint8),  # No deletions
            mask=torch.zeros(dims, dtype=torch.bool),
            sequence_source=torch.full(
                dims,
                fill_value=msa_dataset_source_to_int[MSADataSource.NONE],
                dtype=torch.uint8,
            ),
        )

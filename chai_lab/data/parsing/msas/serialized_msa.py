# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
import dataclasses
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import antipickle
import numpy as np
import pandas as pd

from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.utils.typing import Int32, UInt8, typecheck

logger = logging.getLogger(__name__)


MAPPED_TOKEN_SKIP: Final[int] = -1
MAPPED_TOKEN_INSERTION: Final[int] = -2


@typecheck
@dataclass
class SerializedMSAForSingleSequence:
    """
    Helper class to serialize / deserialize MSAs.

    Serialization uses mix of np.array or polars Series
    (depending on what is more efficient/convenient for every field).
    Polars can contain array fields, but serializes those inefficiently.

    Antipickle serializes both type (+string) and gzips on top.
    """

    query_sequence: str

    species_id: Int32[np.ndarray, "seq"]
    aligned_tokens: UInt8[np.ndarray, "seq tok"]
    deletions: UInt8[np.ndarray, "seq tok"]

    description: pd.Series  # shape: [seq]
    source_database: pd.Series  # shape: [seq]

    def __post_init__(self):
        assert len(self.description) == len(self.source_database) == self.depth
        assert (self.source_database == MSADataSource.QUERY.value).sum() == 1

    @property
    def depth(self) -> int:
        depth, _ = self.aligned_tokens.shape
        return depth

    @property
    def num_tokens(self) -> int:
        _, num_tokens = self.aligned_tokens.shape
        return num_tokens

    def save_to(self, path: Path):
        assert path.suffix == ".gz", "better have it compressed"
        antipickle.dump(dataclasses.asdict(self), filename=path)

    @typecheck
    @classmethod
    def from_query_only(
        cls, query_str: str, query_tokenized: UInt8[np.ndarray, "seq"]
    ) -> "SerializedMSAForSingleSequence":
        query_tokenized = query_tokenized[None, :]  # Insert a seq dim
        deletions = np.zeros_like(query_tokenized)
        species_id = np.zeros((1,), dtype=np.int32)
        description = pd.Series(["query sequence for empty MSA"])
        source_database = pd.Series([MSADataSource.QUERY.value])

        return cls(
            query_sequence=query_str,
            species_id=species_id,
            aligned_tokens=query_tokenized,
            deletions=deletions,
            description=description,
            source_database=source_database,
        )

    @classmethod
    def concatenate(
        cls, serialized_msas: list["SerializedMSAForSingleSequence"]
    ) -> "SerializedMSAForSingleSequence":
        """Concatenate the serialized MSAs along depth dimension."""
        queries = {msa.query_sequence for msa in serialized_msas}
        assert len(queries) == 1

        return cls(
            query_sequence=queries.pop(),
            species_id=np.concatenate([m.species_id for m in serialized_msas], axis=0),
            aligned_tokens=np.concatenate(
                [m.aligned_tokens for m in serialized_msas], axis=0
            ),
            deletions=np.concatenate([m.deletions for m in serialized_msas], axis=0),
            description=pd.concat([m.description for m in serialized_msas]),
            source_database=pd.concat([m.source_database for m in serialized_msas]),
        )

# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
"""
Parsing code for .aligned.pqt files for MSAs.

NOTE this code is meant to be open-sourced, so it should contain as little proprietary
information as possible
"""

import hashlib
from pathlib import Path
from typing import Literal, Mapping

import pandas as pd
import pandera as pa
import torch
from einops import repeat
from pandera.typing import Series

from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.parsing.fasta import read_fasta
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.msas.data_source import (
    MSADataSource,
    msa_dataset_source_to_int,
    msa_dataset_source_to_quota,
)
from chai_lab.data.parsing.msas.species import UNKNOWN_SPECIES, get_tax_ids, stable_hash
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order

RECOGNIZED_SOURCES = {s.value for s in MSADataSource.get_default_sources()}


class AlignedParquetModel(pa.DataFrameModel):
    """Model for aligned parquet files."""

    sequence: Series[str]
    source_database: Series[str] = pa.Field(isin=RECOGNIZED_SOURCES.union({"query"}))
    pairing_key: Series[str]
    comment: Series[str]


def hash_sequence(seq: str) -> str:
    hash_object = hashlib.sha256(seq.encode())
    return hash_object.hexdigest()


def expected_basename(query_sequence: str) -> str:
    """Get the expected filename based on the uppercased query sequence."""
    seqhash = hash_sequence(query_sequence.upper())
    return f"{seqhash}.aligned.pqt"


def _parse_single_source_pqt(
    table: pd.DataFrame, query_seq: str, quota: int | None = None
) -> MSAContext:
    """Parse MSAs for a single source; drops duplicates and adds query as first row."""
    assert table["source_database"].nunique() == 1
    source_parsed = MSADataSource(table["source_database"].values[0])

    if quota is not None:
        assert quota > 0
        table = table.head(quota)

    # Drop duplicates (done after cropping to quote)
    table.drop_duplicates(subset=["sequence"], inplace=True, keep="first")

    # Tokenize the query and the match sequences
    aligned_seqs: list[str] = [query_seq] + table["sequence"].tolist()

    result_sequences, result_deletions = tokenize_sequences_to_arrays(aligned_seqs)
    assert result_sequences.min() >= 0
    assert result_sequences.max() < len(residue_types_with_nucleotides_order)
    depth, n_tokens = result_sequences.shape

    # Build pairing key; first entry corresponds to query which we assign UNKNOWN
    # Empty pairing key is also encoded as UNKNOWN
    tax_id_list = [UNKNOWN_SPECIES] + [
        stable_hash(k) if k else UNKNOWN_SPECIES for k in table["pairing_key"]
    ]
    species = torch.asarray(tax_id_list, dtype=torch.int32)

    tokens = torch.from_numpy(result_sequences)
    species = repeat(species, "seq -> seq n_tokens", n_tokens=n_tokens)
    deletion_matrix = torch.from_numpy(result_deletions)
    mask = torch.ones_like(tokens, dtype=torch.bool)

    return MSAContext(
        dataset_source=source_parsed,
        tokens=tokens,
        species=species,
        deletion_matrix=deletion_matrix,
        mask=mask,
        sequence_source=torch.full_like(
            tokens, fill_value=msa_dataset_source_to_int[source_parsed]
        ),
        is_paired_mask=torch.zeros(depth, dtype=torch.bool),
    )


def parse_aligned_pqt_to_msa_set(
    aligned_pqt_path: Path, apply_quota: bool = True
) -> dict[MSADataSource, MSAContext]:
    """
    Parse .aligned.pqt files, following the schema defined in AlignedParquetModel.
    If apply_quota is specified, then apply per-source quota limits.
    """
    raw_table = pd.read_parquet(aligned_pqt_path)
    try:
        # inplace coerces types in place
        AlignedParquetModel.validate(raw_table, inplace=True)
    except pa.errors.SchemaError as e:
        raise ValueError(f"Invalid schema: {e}")

    # Exactly one query row
    query_row = raw_table.iloc[0]
    assert query_row["source_database"] == "query", "First row must be query"
    query_seq: str = query_row["sequence"]

    # Non-query items
    table = raw_table.iloc[1:]
    retval = {}
    for src, group in table.groupby(by="source_database"):
        assert src != "query", "Encountered query as source"
        parsed_src = MSADataSource(src)
        retval[parsed_src] = _parse_single_source_pqt(
            group,
            query_seq=query_seq,
            quota=msa_dataset_source_to_quota[parsed_src] if apply_quota else None,
        )
    return retval


def a3m_to_aligned_dataframe(
    a3m_path: Path | str,
    source_database: MSADataSource,
    insert_pairing_key: bool = True,
) -> pd.DataFrame:
    """Reformat the a3m as a parquet. Enforces that the first sequence is query."""
    alignments = read_fasta(a3m_path)

    records: list[dict[str, str]] = []
    for i, alignment in enumerate(alignments):
        # Assume first entry in a3m is the query
        src = "query" if i == 0 else source_database.value
        record = {
            "sequence": alignment.sequence,
            "source_database": src,
            # Empty pairing keys get encoded as UNKNOWN
            "pairing_key": str(
                get_tax_ids([alignment.header], source_database)[0]
                if insert_pairing_key
                else ""
            ),
            "comment": "",
        }
        records.append(record)
    assert records[0]["source_database"] == "query"
    retval = pd.DataFrame.from_records(records)
    AlignedParquetModel.validate(retval)
    return retval


def merge_multi_a3m_to_aligned_dataframe(
    msa_a3m_files: Mapping[MSADataSource, Path | str],
    insert_keys_for_sources: Literal["all", "none", "uniprot"] = "uniprot",
) -> pd.DataFrame:
    """Merge multiple a3m files into a single aligned parquet file."""
    dfs = {
        src: a3m_to_aligned_dataframe(
            a3m_path,
            src,
            insert_pairing_key=(
                src == MSADataSource.UNIPROT
                if insert_keys_for_sources == "uniprot"
                else (insert_keys_for_sources == "all")
            ),
        )
        for src, a3m_path in msa_a3m_files.items()
    }
    # Check that all the dfs share the same query sequence
    queries = {df.iloc[0]["sequence"] for df in dfs.values()}
    assert len(queries) == 1
    # As a base, set the query sequence
    chunks = [next(iter(dfs.values())).iloc[0:1]]
    for df in dfs.values():
        # Take the non-query sequences for all sources
        chunks.append(df.iloc[1:])
    return pd.concat(chunks, ignore_index=True).reset_index(drop=True)

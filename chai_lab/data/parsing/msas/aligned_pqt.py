# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
"""
Parsing code for .aligned.pqt files for MSAs.
"""

import hashlib
import logging
from functools import lru_cache
from pathlib import Path
from typing import Literal, Mapping, Optional

import pandas as pd
import pandera as pa
import torch
from einops import repeat
from pandera.typing import Series

from chai_lab.data.dataset.msas.msa_context import NO_PAIRING_KEY, MSAContext
from chai_lab.data.parsing.fasta import read_fasta
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.msas.data_source import (
    MSADataSource,
    msa_dataset_source_to_int,
    msa_dataset_source_to_priority,
    msa_dataset_source_to_quota,
)
from chai_lab.data.parsing.msas.species import get_tax_names
from chai_lab.utils.typing import typecheck

RECOGNIZED_SOURCES: set[str] = {
    s.value for s in MSADataSource.get_default_sources() + [MSADataSource.QUERY]
}


class AlignedParquetModel(pa.DataFrameModel):
    """Model for aligned parquet files."""

    sequence: Series[str]
    source_database: Series[str] = pa.Field(isin=RECOGNIZED_SOURCES)
    pairing_key: Series[str]
    comment: Series[str]


def hash_sequence(seq: str) -> str:
    hash_object = hashlib.sha256(seq.encode())
    return hash_object.hexdigest()


@lru_cache(maxsize=1_000_000)
def stable_hash_for_pairkey(str) -> int:
    # very basic, fast hash, converts to int
    return int(hashlib.sha256(str.encode("utf-8")).hexdigest()[:7], 16)


def expected_basename(query_sequence: str) -> str:
    """Get the expected filename based on the uppercased query sequence."""
    seqhash = hash_sequence(query_sequence.upper())
    return f"{seqhash}.aligned.pqt"


def parse_aligned_pqt_to_msa_context(
    aligned_pqt_path: Path | str,
    quota_sizes: dict[MSADataSource, int] | None = msa_dataset_source_to_quota,
) -> MSAContext:
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

    # Apply quota
    if quota_sizes is not None:
        # For each group, take the quota size if there is a quota found, otherwise take
        # the whole table (some big number)
        # NOTE quota here requires an offset of -1 because in the a3m parsing logic, the
        # quota caps the number of records from a fasta a3m file, which contains the
        # query. Here, the groupby operation does NOT include the query row so we use
        # the offset to account for this.
        raw_table = (
            raw_table.groupby("source_database")
            .apply(
                lambda table: table.head(
                    quota_sizes.get(MSADataSource(table.name), 1_000_000) - 1
                )  # type: ignore
            )
            .reset_index(drop=True)
        )
        assert isinstance(raw_table, pd.DataFrame)

    # Sort by priority, query goes first, other DBs have some order
    raw_table = raw_table.sort_values(
        "source_database",
        key=lambda series: (
            pd.Series(
                [msa_dataset_source_to_priority[MSADataSource(x)] for x in series]
            )
        ),
    )

    sequences, deletions = tokenize_sequences_to_arrays(raw_table["sequence"].tolist())
    _n_seq, n_tok = sequences.shape

    pairing_key_hash = torch.asarray(
        [
            stable_hash_for_pairkey(k) if k else NO_PAIRING_KEY
            for k in raw_table["pairing_key"]
        ],
        dtype=torch.int32,
    )

    source_db_1d = [
        msa_dataset_source_to_int[MSADataSource(x)]
        for x in raw_table["source_database"]
    ]
    source_db_int = repeat(
        torch.asarray(source_db_1d, dtype=torch.uint8), "s -> s t", t=n_tok
    )
    pairing_key_hash = repeat(pairing_key_hash, "s -> s t", t=n_tok)

    tokens = torch.from_numpy(sequences)

    return MSAContext(
        tokens=tokens,
        pairing_key_hash=pairing_key_hash,
        deletion_matrix=torch.from_numpy(deletions),
        sequence_source=source_db_int,
        mask=torch.ones_like(tokens, dtype=torch.bool),
    )


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
        src = MSADataSource.QUERY.value if i == 0 else source_database.value
        record = {
            "sequence": alignment.sequence,
            "source_database": src,
            # Empty pairing keys get encoded as UNKNOWN
            "pairing_key": (
                get_tax_names([alignment.header], source_database).pop()
                if insert_pairing_key
                else ""
            ),
            "comment": alignment.header,
        }
        records.append(record)
    assert records[0]["source_database"] == "query"
    retval = pd.DataFrame.from_records(records)
    AlignedParquetModel.validate(retval)
    return retval


@typecheck
def merge_multi_a3m_to_aligned_dataframe(
    msa_a3m_files: Mapping[Path, MSADataSource],
    insert_keys_for_sources: Literal["all", "none", "uniprot"] = "uniprot",
) -> pd.DataFrame:
    """Merge multiple a3ms from the same query sequence into a single aligned parquet."""
    dfs = {
        src: a3m_to_aligned_dataframe(
            a3m_path,
            src,
            insert_pairing_key=(
                src in (MSADataSource.UNIPROT, MSADataSource.UNIPROT_N3)
                if insert_keys_for_sources == "uniprot"
                else (insert_keys_for_sources == "all")
            ),
        )
        for a3m_path, src in msa_a3m_files.items()
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


def merge_a3m_in_directory(directory: str, output_directory: Optional[str] = None):
    """Finds .a3m files in a directory and combine them into a single aligned.pqt file.
    Files are expected to be named like hits_uniref90.a3m (uniref90 is the source database).
    All files in the directory are assumed to be derived from the same query sequence.

    Provided as a example commandline interface to merge files.
    """
    dir_path = Path(directory)
    assert dir_path.is_dir()

    mapped_a3m_files = {}
    for file in dir_path.glob("*.a3m"):
        # Automatically determine the source based on filename
        dbname = file.stem.replace("_hits", "").replace("hits_", "")
        try:
            msa_src = MSADataSource(dbname)
            logging.info(f"Found {msa_src} MSAs in {file}")
        except ValueError:
            msa_src = MSADataSource.UNIREF90
            logging.warning(
                f"Could not determine source for {file=}; default to {msa_src}"
            )
        mapped_a3m_files[file] = msa_src
    df = merge_multi_a3m_to_aligned_dataframe(
        mapped_a3m_files, insert_keys_for_sources="uniprot"
    )
    # Get the query sequence and use it to determine where we save the file.
    query_seq: str = df.iloc[0]["sequence"]
    # Default to writing into the same directory if output directory isn't specified
    outdir = Path(output_directory) if output_directory is not None else dir_path
    outdir.mkdir(exist_ok=True, parents=True)
    df.to_parquet(outdir / expected_basename(query_seq))


if __name__ == "__main__":
    import typer

    logging.basicConfig(level=logging.INFO)

    typer.run(merge_a3m_in_directory)

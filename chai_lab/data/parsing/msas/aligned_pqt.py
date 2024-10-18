# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
"""
Parsing code for .aligned.pqt files for MSAs.
"""

import hashlib
import logging
from pathlib import Path
from typing import Literal, Mapping

import numpy as np
import pandas as pd
import pandera as pa
from pandera.typing import Series

from chai_lab.data.parsing.fasta import read_fasta
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.msas.data_source import (
    MSADataSource,
    msa_dataset_source_to_priority,
    msa_dataset_source_to_quota,
)
from chai_lab.data.parsing.msas.serialized_msa import SerializedMSAForSingleSequence
from chai_lab.data.parsing.msas.species import (
    UNKNOWN_SPECIES,
    get_tax_names,
    stable_hash,
)
from chai_lab.utils.typing import typecheck

RECOGNIZED_SOURCES: set[str] = {
    s.value for s in MSADataSource.get_default_sources() + [MSADataSource.QUERY]
}


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


def parse_aligned_pqt_to_msa_set(
    aligned_pqt_path: Path | str,
    apply_quota: bool = True,
    quota_sizes: dict[MSADataSource, int] = msa_dataset_source_to_quota,
) -> SerializedMSAForSingleSequence:
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

    # Apply quota
    if apply_quota:
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

    # Sort by priority
    # db2priority is a dict that maps MSADataSource to a non-negative integer priority
    # insert the query as -1 to make sure it is sorted first.
    raw_table = raw_table.sort_values(
        "source_database",
        key=lambda series: (
            pd.Series(
                [
                    (
                        -1
                        if (src := MSADataSource(x)) == MSADataSource.QUERY
                        else msa_dataset_source_to_priority[src]
                    )
                    for x in series
                ]
            )
        ),
    )

    sequences, deletions = tokenize_sequences_to_arrays(raw_table["sequence"].tolist())
    species = np.array(
        [stable_hash(k) if k else UNKNOWN_SPECIES for k in raw_table["pairing_key"]],
        dtype=np.int32,
    )
    return SerializedMSAForSingleSequence(
        query_seq,
        species_id=species,
        aligned_tokens=sequences,
        deletions=deletions,
        description=raw_table["comment"],
        source_database=raw_table["source_database"],
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


def _merge_files_in_directory(directory: str):
    """Finds .a3m files in a directory and combine them into a single aligned.pqt file.
    Files are expected to be named like hits_uniref90.a3m (uniref90 is the source database).
    All files in the directoroy are assumed to be derived from the same query sequence.

    Provided as a example commandline interface to merge files.
    """
    dir_path = Path(directory)
    assert dir_path.is_dir()

    mapped_a3m_files = {}
    for file in dir_path.glob("*.a3m"):
        # Automatically determine the source based on filename
        dbname = file.stem.replace("_hits", "")
        try:
            msa_src = MSADataSource(dbname)
        except Exception:
            msa_src = MSADataSource.UNIREF90
            logging.warning(
                f"Could not determine source for {file=}; default to {msa_src}"
            )
        mapped_a3m_files[msa_src] = file
    df = merge_multi_a3m_to_aligned_dataframe(
        mapped_a3m_files, insert_keys_for_sources="uniprot"
    )
    # Get the query sequence and use it to determine where we save the file.
    query_seq: str = df.iloc[0]["sequence"]
    df.to_parquet(dir_path / expected_basename(query_seq))


if __name__ == "__main__":
    import typer

    typer.run(_merge_files_in_directory)

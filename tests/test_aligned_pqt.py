# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from pathlib import Path

import pandas as pd
import pytest

from chai_lab.data.parsing.msas.aligned_pqt import (
    expected_basename,
    merge_a3m_in_directory,
    merge_multi_a3m_to_aligned_dataframe,
)
from chai_lab.data.parsing.msas.data_source import MSADataSource


@pytest.mark.parametrize(
    "second_source", [MSADataSource.UNIREF90, MSADataSource.MGNIFY]
)
def test_merge_retains_hits_from_every_file(
    tmp_path: Path, second_source: MSADataSource
):
    first = tmp_path / "first.a3m"
    first.write_text(">query\nAG\n>first hit\nAA\n")
    second = tmp_path / "second.a3m"
    second.write_text(">query\nAG\n>second hit\nGG\n")

    table = merge_multi_a3m_to_aligned_dataframe(
        {first: MSADataSource.UNIREF90, second: second_source}
    )

    assert table["sequence"].tolist() == ["AG", "AA", "GG"]
    assert table["source_database"].tolist() == [
        "query",
        MSADataSource.UNIREF90.value,
        second_source.value,
    ]
    assert table["comment"].tolist() == ["query", "first hit", "second hit"]


def test_merge_checks_queries_in_files_from_same_source(tmp_path: Path):
    first = tmp_path / "first.a3m"
    first.write_text(">query\nAG\n>first hit\nAA\n")
    second = tmp_path / "second.a3m"
    second.write_text(">query\nGA\n>second hit\nGG\n")

    with pytest.raises(AssertionError):
        merge_multi_a3m_to_aligned_dataframe(
            {first: MSADataSource.UNIREF90, second: MSADataSource.UNIREF90}
        )


def test_directory_conversion_retains_files_from_same_source(tmp_path: Path):
    (tmp_path / "hits_uniref90.a3m").write_text(">query\nAG\n>first hit\nAA\n")
    (tmp_path / "uniref90_hits.a3m").write_text(">query\nAG\n>second hit\nGG\n")

    merge_a3m_in_directory(str(tmp_path))

    table = pd.read_parquet(tmp_path / expected_basename("AG"))
    assert table["sequence"].iloc[0] == "AG"
    assert set(table["sequence"]) == {"AG", "AA", "GG"}
    assert table["source_database"].tolist().count("query") == 1

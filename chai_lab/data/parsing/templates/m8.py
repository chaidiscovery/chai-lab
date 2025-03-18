# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Iterator

import gemmi
import pandas as pd
import torch

from chai_lab.data.io import rcsb
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.templates.template_hit import TemplateHit
from chai_lab.tools.kalign import kalign_query_to_reference

logger = logging.getLogger(__name__)


def parse_m8_file(fname: Path) -> pd.DataFrame:
    """Parse the m8 alignment format describing template information."""
    table = pd.read_csv(
        fname,
        delimiter="\t",
        header=None,
        names=[
            "query_id",
            "subject_id",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "evalue",
            "bitscore",
            "comment",
        ],
    ).sort_values(by=["query_id", "evalue"])
    for k in ["query_start", "query_end", "subject_start", "subject_end"]:
        table[k] = table[k].astype(int)
    return table


def parse_m8_to_template_hits(
    query_pdb_id: str,
    query_sequence: str,
    m8_path: Path,
    template_cif_folder: Path | None = None,
) -> Iterator[TemplateHit]:
    assert m8_path.is_file()

    if template_cif_folder is not None:
        template_cif_folder.mkdir(parents=True, exist_ok=True)

    table = parse_m8_file(m8_path)

    # Subset to those matching the query pdb id
    table = table.loc[table.query_id.astype(str) == query_pdb_id]
    if len(table) == 0:
        logger.warning(f"[{query_pdb_id=}] No corresponding entries in {m8_path=}")
        return  # No hits

    assert table.query_id.nunique() == 1

    # start and ends are 1-indexed in the input table
    # update the starts to be 0 indexed w.t. the indexing becomes 0-index half open
    for k in ["query_start", "subject_start"]:
        table[k] = table[k] - 1
        assert (table[k] >= 0).all()

    # Process each hit
    counter = 0
    for i, row in enumerate(table.itertuples()):
        hit_identifier, hit_chain = row.subject_id.split("_")  # type: ignore
        assert isinstance(hit_identifier, str) and isinstance(hit_chain, str)
        with TemporaryDirectory() as tmpdir:
            cif_file = rcsb.download_cif_file(
                hit_identifier.upper(),
                directory=(
                    Path(tmpdir) if template_cif_folder is None else template_cif_folder
                ),
            )
            structure = gemmi.read_structure(path=str(cif_file))

        chain: gemmi.Chain = structure[0][hit_chain]  # Indexes by auth chain
        # NOTE this sequence excludes unresolved residues, and adds "-" to indicate gap
        # We strip the gap because the m8 file is indexed without the gap
        seq: str = chain.get_polymer().make_one_letter_sequence().replace("-", "")

        seq_matching = seq[row.subject_start : row.subject_end]  # type: ignore

        alignment = kalign_query_to_reference(ref=query_sequence, query=seq_matching)
        if alignment is None:
            continue

        match_seq, match_del = tokenize_sequences_to_arrays([alignment.query_a3m_line])

        try:
            hit = TemplateHit(
                query_pdb_id=query_pdb_id,
                query_sequence=query_sequence,
                index=i,
                pdb_id=hit_identifier,
                chain_id=hit_chain,
                hit_start=row.subject_start,  # type: ignore
                hit_end=row.subject_end,  # type: ignore
                hit_tokens=torch.from_numpy(match_seq[0]).int(),
                deletion_matrix=torch.from_numpy(match_del[0]),
                query_seq_realigned=alignment.query_a3m_line,
                cif_path=cif_file if template_cif_folder is not None else None,
            )
            assert hit.query_start_end == alignment.reference_span

            logger.info(f"[{query_pdb_id}] {hit}")
            yield hit

            counter += 1

        except Exception:
            logger.warning(
                f"[{query_pdb_id=}] Could not load template from {hit_identifier} {hit_chain}",
                exc_info=True,
            )
            pass

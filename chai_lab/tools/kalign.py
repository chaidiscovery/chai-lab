# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.


import logging
import shutil
import subprocess
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory

from chai_lab.data.parsing.fasta import Fasta, read_fasta, write_fastas

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class KalignAlignment:
    """Alignment expressing how a given query matches to a reference."""

    reference_aligned: str  # the aligned output of kalign
    query_aligned: str  # aligned output of kalign

    def __post_init__(self):
        assert len(self.reference_aligned) == len(self.query_aligned) > 0
        # There should be no positions where both the query and reference are gaps
        assert not any(
            r == "-" and q == "-"
            for r, q in zip(self.reference_aligned, self.query_aligned)
        )

    @property
    def query_a3m_line(self) -> str:
        # If the reference is a "-" that indicates that the corresponding query residue
        # is an insertion; mark it as such
        retval = [
            query_res if ref_res != "-" else query_res.lower()
            for query_res, ref_res in zip(self.query_aligned, self.reference_aligned)
        ]
        return "".join(retval)

    @property
    def reference_span(self) -> tuple[int, int]:
        """The span of the reference that is matched by the query.

        Includes middle gaps, but excludes trailing gaps."""
        reference_positions_covered = set()
        i = 0  # Reference pointer
        for query_res, ref_res in zip(self.query_aligned, self.reference_aligned):
            if query_res != "-":  # Query is not a gap; matched
                reference_positions_covered.add(i)
            if ref_res != "-":  # Move pointer along if the reference is not gapped
                i += 1
        reference_span = sorted(reference_positions_covered)
        return reference_span[0], reference_span[-1]


# NOTE caching here helps in the case of dimers that have the same chain sequence that
# we process twice in a row. This allows us to not re-compute alignments.
@lru_cache(maxsize=1024)
def kalign_query_to_reference(
    ref: str, query: str, threads: int = 1
) -> KalignAlignment | None:
    """Align the query to the reference; ignores all insertions (lower case) and
    deletions (- gaps) in the query. This happens by casting all inputs to upper case
    and removing all - chars.

    Returns KalignAlignment object that encapsulates the aligned query w.r.t. reference
    """
    assert (
        shutil.which("kalign") is not None
    ), "kalign is required for templates, but was not found in PATH. Try 'apt install kalign'?"

    query = query.upper().replace("-", "")
    with TemporaryDirectory() as tmp_dir:
        # Write the input sequences
        fastas = [
            Fasta("ref", ref),
            Fasta("query", query),
        ]
        fasta_path = Path(tmp_dir) / "kalign_input.fasta"
        write_fastas(fastas, output_path=str(fasta_path))
        out_path = Path(tmp_dir) / "kalign_output.fasta"
        command = [
            "kalign",
            "-i",
            str(fasta_path),
            "-o",
            str(out_path),
            "--nthreads",
            str(threads),
        ]
        logger.debug(f"Running kalign: {command}")
        try:
            subprocess.run(
                command,
                check=True,
                close_fds=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            logger.debug("Kalign called successfully")
        except subprocess.CalledProcessError:
            logger.error("Error running kalign")
            return None

        aligned_ref, aligned_query = read_fasta(out_path)
        # Check that the output is the same order as input
        assert aligned_ref.header == "ref" and aligned_query.header == "query"

    return KalignAlignment(
        reference_aligned=aligned_ref.sequence, query_aligned=aligned_query.sequence
    )

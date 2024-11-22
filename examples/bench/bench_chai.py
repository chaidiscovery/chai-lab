# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

"""
Helper functions to bench chai1
"""

import logging
from pathlib import Path

import pandas as pd
import torch

from chai_lab import chai1
from chai_lab.data.parsing.fasta import read_fasta
from chai_lab.metrics.dockq import calc_dockq


def run_chai1(input_file: Path | str, outdir: Path) -> dict[Path, float]:
    input_file = Path(input_file)
    outdir.mkdir(exist_ok=True, parents=True)
    structure_candidates = chai1.run_inference(
        fasta_file=input_file,
        output_dir=outdir / input_file.stem,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=42,
        device=torch.device("cuda:0"),
        use_esm_embeddings=True,
    ).sort_by_rank()
    return {
        p: s
        for p, s in zip(
            structure_candidates.cif_paths,
            structure_candidates.candidate_aggregate_scores,
        )
    }


def bench_multimers(input_base_dir: Path, reference_base_dir: Path) -> pd.DataFrame:
    inputs: list[Path] = []
    for fasta in input_base_dir.glob("*.fasta"):
        protein_chains = [
            chain for chain in read_fasta(fasta) if chain.header.startswith("protein")
        ]
        if len(protein_chains) < 2:
            continue
        inputs.append(fasta)

    # For each input, run chai1
    records = []
    for fasta in inputs:
        # Determine the reference file
        ref_paths = (
            list(reference_base_dir.glob(f"{fasta.stem}*.cif"))
            + list(reference_base_dir.glob(f"{fasta.stem}*.pdb"))
            + list(reference_base_dir.glob(f"{fasta.stem}*.cif.gz"))
        )
        if len(ref_paths) != 1:
            logging.warning(f"[{fasta}] No unique reference structure, skipping")
            continue
        ref_path = ref_paths.pop()

        chai_results = run_chai1(
            fasta, outdir=Path(__file__).parent / "outputs/chai" / fasta.stem
        )
        for rank, pred_path in enumerate(chai_results.keys()):
            dockq_results = calc_dockq(model=pred_path, native=ref_path)
            print(dockq_results)
            records.append({"rank": rank, "dockQ": dockq_results["best_dockq"]})
    return pd.DataFrame.from_records(records)


if __name__ == "__main__":
    basedir = Path(__file__).parent
    df = bench_multimers(
        basedir / "data/casp15/inputs/",
        basedir / "data/casp15/references/",
    )
    print(df)

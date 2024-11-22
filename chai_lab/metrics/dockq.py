# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

import itertools
import json
import logging
from functools import partial
from typing import Any

import typer
from DockQ.DockQ import (
    count_chain_combinations,
    format_mapping,
    format_mapping_string,
    get_all_chain_maps,
    group_chains,
    load_PDB,
    run_on_all_native_interfaces,
)
from parallelbar import progress_map


def calc_dockq(
    model: str,
    native: str,
    mapping: str = "",
    capri_peptide: bool = False,
    small_molecule: bool = False,
    n_cpu: int = 8,
    max_chunk: int = 512,
    allowed_mismatches: int = 0,
    json_out: str = "",
) -> dict[str, Any] | None:
    """
    Lightly modified from the main function in the official DockQ implementation.
    """

    initial_mapping, model_chains, native_chains = format_mapping(
        mapping, small_molecule
    )
    model_structure = load_PDB(
        model, chains=model_chains, small_molecule=small_molecule
    )
    native_structure = load_PDB(
        native, chains=native_chains, small_molecule=small_molecule
    )
    # check user-given chains are in the structures
    model_chains = [c.id for c in model_structure] if not model_chains else model_chains
    native_chains = (
        [c.id for c in native_structure] if not native_chains else native_chains
    )

    if len(model_chains) < 2 or len(native_chains) < 2:
        logging.warning("Need at least two chains in the two inputs")
        return None

    # permute chains and run on a for loop
    best_dockq = -1
    best_result = None
    best_mapping = None

    model_chains_to_combo = [
        mc for mc in model_chains if mc not in initial_mapping.values()
    ]
    native_chains_to_combo = [
        nc for nc in native_chains if nc not in initial_mapping.keys()
    ]

    chain_clusters, reverse_map = group_chains(
        model_structure,
        native_structure,
        model_chains_to_combo,
        native_chains_to_combo,
        allowed_mismatches,
    )
    chain_maps = get_all_chain_maps(
        chain_clusters,
        initial_mapping,
        reverse_map,
        model_chains_to_combo,
        native_chains_to_combo,
    )

    num_chain_combinations = count_chain_combinations(chain_clusters)
    # copy iterator to use later
    chain_maps, chain_maps_ = itertools.tee(chain_maps)

    low_memory = num_chain_combinations > 100
    run_chain_map = partial(
        run_on_all_native_interfaces,
        model_structure,
        native_structure,
        no_align=False,
        capri_peptide=capri_peptide,
        low_memory=low_memory,
    )

    if num_chain_combinations > 1:
        cpus = min(num_chain_combinations, n_cpu)
        chunk_size = min(max_chunk, max(1, num_chain_combinations // cpus))

        # for large num_chain_combinations it should be possible to divide the chain_maps in chunks
        result_this_mappings = progress_map(
            run_chain_map,
            chain_maps,
            total=num_chain_combinations,
            n_cpu=cpus,
            chunk_size=chunk_size,
        )

        for chain_map, (result_this_mapping, total_dockq) in zip(
            chain_maps_, result_this_mappings
        ):
            if total_dockq > best_dockq:
                best_dockq = total_dockq
                best_result = result_this_mapping
                best_mapping = chain_map

        if low_memory:  # retrieve the full output by rerunning the best chain mapping
            best_result, total_dockq = run_on_all_native_interfaces(
                model_structure,
                native_structure,
                chain_map=best_mapping,
                no_align=False,
                capri_peptide=capri_peptide,
                low_memory=False,
            )

    else:  # skip multi-threading for single jobs (skip the bar basically)
        best_mapping = next(chain_maps)
        best_result, best_dockq = run_chain_map(best_mapping)

    if not best_result:
        logging.error(
            "Could not find interfaces in the native model. Please double check the inputs or select different chains with the --mapping flag."
        )
        sys.exit(1)

    info = dict()
    info["model"] = model
    info["native"] = native
    info["best_dockq"] = best_dockq
    info["best_result"] = best_result
    info["GlobalDockQ"] = best_dockq / len(best_result)
    info["best_mapping"] = best_mapping
    info["best_mapping_str"] = f"{format_mapping_string(best_mapping)}"

    if json_out:
        with open(json_out, "w") as sink:
            json.dump(info, sink, indent=2)

    return info


if __name__ == "__main__":
    typer.run(calc_dockq)

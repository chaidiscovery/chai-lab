# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
"""
Code for loading MSAs given
"""

import logging
from pathlib import Path

import torch

from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.dataset.msas.preprocess import (
    concatenate_paired_and_main_msas,
    drop_duplicates,
    merge_main_msas_by_chain,
    merge_msas_by_datasource,
    pair_msas_by_chain_with_species_matching,
    partition_msa_by_pairing_key,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.msas.aligned_pqt import (
    expected_basename,
    parse_aligned_pqt_to_msa_set,
)
from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.data.parsing.structure.entity_type import EntityType


def get_msa_contexts(
    chains: list[Chain],
    msa_directory: Path,
) -> tuple[MSAContext, MSAContext]:
    """Returns two contexts

    - First context to tokenize and give to model
    - Second context for computing summary statistics.
    """
    pdb_ids = set(chain.entity_data.pdb_id for chain in chains)
    assert len(pdb_ids) == 1, f"Found >1 pdb ids in chains: {pdb_ids=}"

    # MSAs are constructed based on sequence, so use the unique sequences present
    # in input chains to determine the MSAs that need to be loaded
    msas_to_load = {
        (chain.entity_data.sequence, chain.entity_data.entity_id)
        for chain in chains
        if chain.entity_data.entity_type == EntityType.PROTEIN
    }

    # Load up the MSAs for each chain; do this by checking for MSA sequences
    msa_contexts_for_entities = dict()
    for seq, entity_id in msas_to_load:
        path = msa_directory / expected_basename(seq)
        if not path.is_file():
            logging.warning(f"No MSA found for sequence: {seq}")
            continue
        msa_contexts_for_entities[(seq, entity_id)] = parse_aligned_pqt_to_msa_set(path)

    # For each chain, either fetch the corresponding MSA or create an empty MSA if it is missing
    msa_sets = [
        (
            msa_contexts_for_entities[k]
            if (k := (chain.entity_data.sequence, chain.entity_data.entity_id))
            in msa_contexts_for_entities
            else {
                MSADataSource.NONE: MSAContext.create(
                    MSADataSource.NONE,
                    tokens=torch.from_numpy(
                        tokenize_sequences_to_arrays([chain.entity_data.sequence])[
                            0
                        ].squeeze(0)
                    ),
                )
            }
        )
        for chain in chains
    ]

    # Stack the MSA for each chain together across MSA sources
    msa_contexts = [merge_msas_by_datasource(msa_set) for msa_set in msa_sets]

    # Re-index to handle residues that are tokenized per-atom
    msa_sets_exploded = [
        msa_ctx[..., chain.structure_context.token_residue_index]
        for chain, msa_ctx in zip(chains, msa_contexts, strict=True)
    ]

    # Pair up the MSAs that have a pairing key (typically species) provided
    divided = [partition_msa_by_pairing_key(m) for m in msa_sets_exploded]
    pairing_contexts = [d[0] for d in divided]
    paired_msa = pair_msas_by_chain_with_species_matching(pairing_contexts)

    # Process main MSA - deduplicate and merge across chains
    main_contexts = [d[1] for d in divided]
    main_msa_deduped = [drop_duplicates(msa) for msa in main_contexts]
    main_msa = merge_main_msas_by_chain(main_msa_deduped)

    # Combine the paired and main MSAs
    merged_msa = concatenate_paired_and_main_msas(paired_msa, main_msa)
    merged_dedup_msa = drop_duplicates(merged_msa)

    return merged_dedup_msa, main_msa

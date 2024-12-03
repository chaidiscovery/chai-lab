# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from pathlib import Path

import torch

from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.dataset.msas.preprocess import (
    drop_duplicates,
    merge_main_msas_by_chain,
    pair_and_merge_msas,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.msas.a3m import tokenize_sequences_to_arrays
from chai_lab.data.parsing.msas.aligned_pqt import (
    expected_basename,
    parse_aligned_pqt_to_msa_context,
)
from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.data.parsing.structure.entity_type import EntityType

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_msa_contexts(
    chains: list[Chain],
    msa_directory: Path,
) -> tuple[MSAContext, MSAContext]:
    """
    Looks inside msa_directory to find .aligned.pqt files to load alignments from.

    Returns two contexts

    - First context to tokenize and give to model
    - Second context for computing summary statistics
    """

    pdb_ids = set(chain.entity_data.pdb_id for chain in chains)
    assert len(pdb_ids) == 1, f"Found >1 pdb ids in chains: {pdb_ids=}"

    # MSAs are constructed based on sequence, so use the unique sequences present
    # in input chains to determine the MSAs that need to be loaded
    def get_msa_contexts_for_seq(seq: str, etype: EntityType) -> MSAContext:
        path = msa_directory / expected_basename(seq)
        # If the MSA is missing, or the query is not a protein, return an empty MSA
        if not path.is_file() or etype != EntityType.PROTEIN:
            if etype == EntityType.PROTEIN:
                # Warn for proteins that have missing MSAs
                logger.warning(f"No MSA found for sequence: {seq}")
            [tokenized_seq], _ = tokenize_sequences_to_arrays([seq])
            return MSAContext.create_single_seq(
                MSADataSource.QUERY, tokens=torch.from_numpy(tokenized_seq)
            )
        msa = parse_aligned_pqt_to_msa_context(path)
        logger.info(f"MSA found for sequence: {seq}, {msa.depth=}")
        return msa

    # For each chain, either fetch the corresponding MSA or create an empty MSA if it is missing
    # + reindex to handle residues that are tokenized per-atom (this also crops if necessary)
    msa_contexts = [
        get_msa_contexts_for_seq(
            seq=chain.entity_data.sequence, etype=chain.entity_data.entity_type
        )[:, chain.structure_context.token_residue_index]
        for chain in chains
    ]

    # used later only for profile statistics
    profile_msa = merge_main_msas_by_chain(
        [drop_duplicates(msa) for msa in msa_contexts]
    )

    joined_msa = pair_and_merge_msas(msa_contexts)
    joined_msa = drop_duplicates(joined_msa)  # rare dups after pairings

    logger.info(f"Prepared MSA context with {joined_msa.depth=}")
    return joined_msa, profile_msa

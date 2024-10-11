# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
import logging
from collections import defaultdict

import torch
from einops import rearrange, reduce
from torch import Tensor

from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.parsing.msas.data_source import MSADataSource
from chai_lab.data.parsing.msas.species import UNKNOWN_SPECIES
from chai_lab.utils.tensor_utils import unique_indexes
from chai_lab.utils.typing import Bool, Int, typecheck

MAX_PAIRED_DEPTH = 8_191
FULL_DEPTH = 16_384


def merge_msas_by_datasource(
    msas: dict[MSADataSource, MSAContext],
    new_source: MSADataSource = MSADataSource.MAIN,
) -> MSAContext:
    """Merges MSAs by data source, concatenating over depth dimension."""
    filtered_first_msa_rows = [
        msa.tokens[0, :] for msa in msas.values() if msa.depth > 0
    ]
    query_sequences_matches = [
        torch.equal(filtered_first_msa_rows[0], row) for row in filtered_first_msa_rows
    ]

    assert all(query_sequences_matches)

    # Only keep query sequence for the first MSA.
    # If an MSA has depth 1 (just query sequence), or no sequence, it is dropped
    merged_msa = MSAContext.cat(
        [
            msa if idx == 0 else msa[1:, :]
            for idx, msa in enumerate(msas.values())
            if msa.depth > 1
        ],
        dataset_source=new_source,
        dim=0,
    )
    if merged_msa.depth == 0:
        # msas from all datasources are empty, use dummy MSA
        return MSAContext.create_empty(merged_msa.num_tokens)
    return merged_msa


def merge_main_msas_by_chain(msas: list[MSAContext]) -> MSAContext:
    """Merges MSAs across chains, concatenating over token dimension."""
    # Skip over fully cropped MSAs
    msas = [msa for msa in msas if msa.num_tokens > 0]

    # Make sure the MSAs all have the same depth
    max_msa_depth = max(msa.depth for msa in msas)
    padded_msas = [msa.pad(max_msa_depth=max_msa_depth) for msa in msas]

    # Concatenate the MSAs along the tokens dimension
    return MSAContext.cat(
        padded_msas,
        dataset_source=MSADataSource.MAIN,
        dim=1,
    )


def partition_msa_by_pairing_key(msa: MSAContext) -> tuple[MSAContext, MSAContext]:
    """Divide the given msa by whether or not a pairing key is set."""
    pairing_key_set = msa.species != 0
    return msa[pairing_key_set, ...], msa[~pairing_key_set, ...]


@typecheck
def _pair_msas_by_chain_with_species_matching_just_pairing(
    msas: list[MSAContext],  # msas from different chains
) -> tuple[
    Int[Tensor, "num_msas sum_num_indices"], Bool[Tensor, "num_msas sum_num_indices"]
]:
    """Pairs based on species matrix values."""
    # for each species, dict keyed by msa index
    # containing a list of sequence indices that are from the given species
    species2msa_id2sequence: dict[str, dict[int, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for msa_id, msa in enumerate(msas):
        for i, species_vec in enumerate(msa.species):
            # input MSA should only have one species per row
            uniq_species = torch.unique(species_vec)
            assert uniq_species.numel() == 1
            species = uniq_species.item()
            # skip unknown species and skip reference sequence
            if species != UNKNOWN_SPECIES and i != 0:
                species2msa_id2sequence[species][msa_id].append(i)

    # ensures correct shape if no other chunks added
    dummy_chunk = torch.zeros([len(msas), 0], dtype=torch.long)
    chunks = [dummy_chunk]
    masks = [torch.zeros([len(msas), 0], dtype=torch.bool)]

    # Hack: When we pair UniProt and UniProt-n3, we may have different query sequences
    # (we are inconsistent in how we deal with modified amino acids) which therefore
    # results in an MSA with depth 2. Let's skip pairing in these cases.
    num_msas_to_pair = len([msa for msa in msas if msa.depth > 2])

    # For each MSA, precompute similarities between each row and the reference sequence
    msa_similarity = [
        reduce(msa.tokens[0, :] != msa.tokens, "seq tok -> seq", "sum") for msa in msas
    ]

    for _species, msa_id2sequence in species2msa_id2sequence.items():
        if len(msa_id2sequence) < num_msas_to_pair:
            continue  # not present in one of non-empty MSAs

        # TODO current strategy works poorly for complexes with > 2 proteins,
        #  if one chain misses a species, pairing of all others won't be considered
        n_pairings = min(len(seqs) for seqs in msa_id2sequence.values())

        # select closest by edit distance
        chunk = torch.zeros([len(msas), n_pairings], dtype=torch.long)
        mask = torch.zeros([len(msas), n_pairings], dtype=torch.bool)

        for msa_id, msa in enumerate(msas):
            if msa_id in msa_id2sequence:
                seq_idx = torch.asarray(msa_id2sequence[msa_id], dtype=torch.long)
                edit_dists = msa_similarity[msa_id][seq_idx]

                chunk[msa_id, :] = seq_idx[torch.argsort(edit_dists)][:n_pairings]
                mask[msa_id, :] = True

        chunks.append(chunk)
        masks.append(mask)

    # num_msas x sum_num_indices
    return torch.concatenate(chunks, dim=1), torch.concatenate(masks, dim=1)


def pair_msas_by_chain_with_species_matching(msas: list[MSAContext]) -> MSAContext:
    """
    Combines the input MSAs (from different chains) across the different chains into a
    single MSA, merging chains along the tokens dimension and pairing by pairing key.
    """
    if len(msas) == 0:
        logging.warning("No MSAs to pair; they are all empty")
        return MSAContext.create_empty(0)
    assert all(msa.num_tokens > 0 for msa in msas)

    # Create paired part of the MSA
    all_stacked_matches, all_stacked_masks = [
        x[:, :MAX_PAIRED_DEPTH]
        for x in _pair_msas_by_chain_with_species_matching_just_pairing(msas)
    ]
    num_msas, sum_num_indices = all_stacked_matches.shape
    assert num_msas == len(msas)

    # Truncate the MSAs to the query sequence
    query_msas = [msa[0:1, :] for msa in msas]

    # Get the paired sequences for each MSA (if we found any species matches)
    paired_msas = []
    for msa_idx, (matches, mask) in enumerate(
        zip(all_stacked_matches, all_stacked_masks, strict=True)
    ):
        msa_chunk = msas[msa_idx][matches, :]
        mask = rearrange(mask, "depth -> depth 1")
        paired_msas.append(msa_chunk.apply_mask(mask))

    # Stack the query and paired MSAs along the depth dimension
    full_msas = [
        MSAContext.cat(
            [query_msa, paired_msa],
            dataset_source=MSADataSource.PAIRED,
            dim=0,
        )
        for query_msa, paired_msa in zip(query_msas, paired_msas)
    ]

    # Finally, merge the MSAs into a single MSA along the tokens dimension
    merged_msas = MSAContext.cat(
        full_msas,
        dataset_source=MSADataSource.PAIRED,
        dim=1,
    )
    assert merged_msas.depth == sum_num_indices + 1
    assert merged_msas.num_tokens == sum(msa.num_tokens for msa in msas)

    merged_msas.is_paired_mask = torch.ones_like(
        merged_msas.is_paired_mask,
        dtype=torch.bool,
    )

    return merged_msas


def drop_duplicates(msa: MSAContext) -> MSAContext:
    # when chains are cropped out, they have zero tokens. We should not try to
    # deduplicate in this case.
    if msa.tokens.shape[-1] == 0:
        return msa
    # sorted: don't shuffle the paired and main MSA
    _, idxs = unique_indexes(msa.tokens, dim=0, sorted=True)
    return msa[idxs, :]


@typecheck
def concatenate_paired_and_main_msas(
    paired_msa: MSAContext,
    main_msa: MSAContext,
) -> MSAContext:
    """
    The paired MSA and main MSA can have different depths.
    """
    paired_depth = min(paired_msa.depth, MAX_PAIRED_DEPTH)
    truncated_paired_msa = paired_msa[:paired_depth, :]

    main_depth = min(main_msa.depth, FULL_DEPTH - paired_depth)
    truncated_main_msa = main_msa[:main_depth, :]

    return MSAContext.cat(
        [truncated_paired_msa, truncated_main_msa],
        dim=0,
    )

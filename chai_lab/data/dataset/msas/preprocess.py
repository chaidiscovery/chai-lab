# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
import logging
from collections import Counter, defaultdict

import torch
from einops import reduce
from torch import Tensor

from chai_lab.data.dataset.msas.msa_context import NO_PAIRING_KEY, MSAContext
from chai_lab.utils.tensor_utils import unique_indexes
from chai_lab.utils.typing import Int, typecheck

MAX_PAIRED_DEPTH = 8_192
FULL_DEPTH = 16_384
_UKEY_FOR_QUERY = (-999, -999)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def merge_main_msas_by_chain(msas: list[MSAContext]) -> MSAContext:
    """Merges MSAs across chains, concatenating over token dimension."""
    # Skip over fully cropped MSAs
    msas = [msa for msa in msas if msa.num_tokens > 0]

    # Make sure the MSAs all have the same depth
    max_msa_depth = max(msa.depth for msa in msas)
    padded_msas = [msa.pad(max_msa_depth=max_msa_depth) for msa in msas]

    # Concatenate the MSAs along the tokens dimension
    return MSAContext.cat(padded_msas, dim=1)


def drop_duplicates(msa: MSAContext) -> MSAContext:
    # when chains are cropped out, they have zero tokens. We should not try to
    # deduplicate in this case.
    if msa.tokens.shape[-1] == 0:
        return msa
    if msa.depth == 0:
        return msa
    # sorted: don't shuffle the paired and main MSA
    _, idxs = unique_indexes(msa.tokens, dim=0, sorted=True)
    return msa[idxs, :]


@typecheck
def prepair_ukey(
    pairing_keys: Int[Tensor, "depth"], edit_distances: Int[Tensor, "depth"]
) -> dict[tuple[int, int], int]:
    """Return mapping of (pair_key, rank) -> ith row"""
    # pairing key is not unique, so we handle this with ukey,
    # which is unique and would be mapped 1-to-1.
    # e.g. keys 1, 1, 2, 1, 2, 3 would be transformed into tuples:
    # their ukeys = (1, 1), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1)

    # For each row, find its rank by edit distance within hits sharing pairing key
    key2edits: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for idx, (pairing_key, edit_distance) in enumerate(
        zip(pairing_keys.tolist(), edit_distances.tolist())
    ):
        if idx == 0:  # Do not count the query
            continue
        key2edits[pairing_key].append((edit_distance, idx))
    idx2pairrank: dict[int, tuple[int, int]] = {}
    for pairkey, row_edits in key2edits.items():
        edits, row_indices = zip(*row_edits)
        assert len(edits) == len(row_indices)
        ranks = torch.argsort(torch.tensor(edits))
        for idx, rank in zip(row_indices, ranks, strict=True):
            idx2pairrank[idx] = (pairkey, rank.item())

    result: dict[tuple[int, int], int] = {}
    for rowid, pairkey in enumerate(pairing_keys.tolist()):
        if rowid == 0:
            result[_UKEY_FOR_QUERY] = rowid  # special ukey for query
            continue
        if pairkey == NO_PAIRING_KEY:
            continue  # should not be paired
        ukey = idx2pairrank[rowid]
        result[ukey] = rowid
    return result


def pair_and_merge_msas(msas: list[MSAContext]) -> MSAContext:
    msa_edit_distances = [
        reduce(msa.tokens[0, :] != msa.tokens, "seq tok -> seq", "sum") for msa in msas
    ]
    chain2ukey2rowid = [
        prepair_ukey(msa.pairing_key_hash[:, 0], msa_edit_distance)
        for msa, msa_edit_distance in zip(msas, msa_edit_distances, strict=True)
    ]

    _ukey2count = Counter(
        ukey for ukey2rowid in chain2ukey2rowid for ukey in ukey2rowid.keys()
    )
    # to be paired, key should be in all chains with non-empty msa
    n_msas = sum(msa.depth > 1 for msa in msas)
    selected_ukeys = [ukey for ukey, count in _ukey2count.items() if count >= n_msas]
    selected_ukeys = selected_ukeys[:MAX_PAIRED_DEPTH]
    logger.info(f"Paired {len(selected_ukeys)} rows")
    assert selected_ukeys[0] == _UKEY_FOR_QUERY

    reordered_msas: list[MSAContext] = []

    for ukey2rowid, msa in zip(chain2ukey2rowid, msas, strict=True):
        paired_rowids = [ukey2rowid.get(ukey, None) for ukey in selected_ukeys]
        assert paired_rowids[0] == 0
        _paired_rowids_set = set(paired_rowids)
        assert set(paired_rowids)

        all_rowids = paired_rowids + [
            i for i in range(msa.depth) if i not in _paired_rowids_set
        ]
        all_rowids = all_rowids[:FULL_DEPTH]

        n_paired_msa = sum(x is not None for x in paired_rowids)
        n_unpaired_msa = sum(x is not None for x in all_rowids) - n_paired_msa

        selected_msa = msa.take_rows_with_padding(all_rowids)

        logger.info(
            f"Loaded (paired includes query sequence): "
            f"{n_paired_msa=} {n_unpaired_msa=} out of {msa.depth=} "
        )

        # reorder each msa to have paired elements first; that's same # of rows
        reordered_msas.append(selected_msa)

    result = merge_main_msas_by_chain(reordered_msas)  # pad and concat

    return result

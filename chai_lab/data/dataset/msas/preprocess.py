# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.
from collections import Counter
from typing import Iterable

from chai_lab.data.dataset.msas.msa_context import NO_PAIRING_KEY, MSAContext
from chai_lab.utils.tensor_utils import unique_indexes

MAX_PAIRED_DEPTH = 8_192
FULL_DEPTH = 16_384


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


def pair_and_merge_msas(msas: list[MSAContext]) -> MSAContext:
    # pairing key is not unique, so we handle this with ukey,
    # which is unique and would be mapped 1-to-1.
    # e.g. keys 1, 1, 2, 1, 2, 3 would be transformed into tuples:
    # their ukeys = (1, 1), (1, 2), (2, 1), (1, 3), (2, 2), (3, 1)
    _UKEY_FOR_QUERY = (-999, -999)

    def prepair_ukey(pairing_keys: Iterable[int]) -> dict[tuple[int, int], int]:
        pairkey2count: Counter[int] = Counter()
        result = {}
        for rowid, pairkey in enumerate(pairing_keys):
            if rowid == 0:
                result[_UKEY_FOR_QUERY] = rowid  # special ukey for query
                continue
            if pairkey == NO_PAIRING_KEY:
                continue  # should not be paired
            pairkey2count[pairkey] += 1
            ukey = (pairkey, pairkey2count[pairkey])
            result[ukey] = rowid

        return result

    chain2ukey2rowid = [
        prepair_ukey(msa.pairing_key_hash[:, 0].tolist()) for msa in msas
    ]

    _ukey2count = Counter(
        ukey for ukey2rowid in chain2ukey2rowid for ukey in ukey2rowid.keys()
    )
    # to be paired, key should be in all chains with non-empty msa
    n_msas = sum(msa.depth > 1 for msa in msas)
    selected_ukeys = [ukey for ukey, count in _ukey2count.items() if count >= n_msas]
    selected_ukeys = selected_ukeys[:MAX_PAIRED_DEPTH]
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

        print(
            f"Loaded (paired in includes query sequence):"
            f"\n{n_paired_msa=} {n_unpaired_msa=} out of {msa.depth=} "
        )

        # reorder each msa to have paired elements first; that's same # of rows
        reordered_msas.append(selected_msa)

    result = merge_main_msas_by_chain(reordered_msas)  # pad and concat

    return result

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
import torch

from chai_lab.data.dataset.msas.msa_context import NO_PAIRING_KEY
from chai_lab.data.dataset.msas.preprocess import _UKEY_FOR_QUERY, prepair_ukey


def test_prepair_ukey():
    keys = torch.tensor([1, 1, 2, 1, NO_PAIRING_KEY, 2, 3])
    edit_dists = torch.arange(len(keys))

    paired = prepair_ukey(keys, edit_dists)
    assert list(paired) == [_UKEY_FOR_QUERY, (1, 0), (2, 0), (1, 1), (2, 1), (3, 0)]
    assert set(paired.values()) == set(
        [i for i, val in enumerate(keys.tolist()) if val != NO_PAIRING_KEY]
    )

    # Reverse the edit distances
    paired = prepair_ukey(keys, torch.tensor(edit_dists.tolist()[::-1]))
    assert list(paired) == [_UKEY_FOR_QUERY, (1, 1), (2, 1), (1, 0), (2, 0), (3, 0)]
    assert set(paired.values()) == set(
        [i for i, val in enumerate(keys.tolist()) if val != NO_PAIRING_KEY]
    )

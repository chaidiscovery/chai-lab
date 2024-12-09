# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

from typing import TypeVar

K = TypeVar("K")
V = TypeVar("V")


def list_dict_to_dict_list(list_dict: list[dict[K, V]]) -> dict[K, list[V]]:
    """
    Converts a list of dicts that contain the same keys to a dict of lists, where each
    list contains an ordered list of values of the corresponding dict.
    """
    if len(list_dict) == 0:
        return {}

    keys = list_dict[0].keys()
    if any(d.keys() != keys for d in list_dict):
        raise ValueError("All dicts must have the same keys")

    return {k: [d[k] for d in list_dict] for k in keys}

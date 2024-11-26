# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Utility classes and functions for feature representations"""

from chai_lab.utils.typing import typecheck


@typecheck
def get_entry_for_key(data: dict, key: str):
    """finds entry 'key' in data dictionary

    Parameters:
        data: the dict to search in
        key: the key to search for

    Example 1:
        data=dict(foo=dict(bar="bar"))
        key = "foo"
        returns: dict(bar="bar")
    Example 2:
        data=dict(foo=dict(bar="bar"))
        key = "foo/bar"
        returns: "bar"

    """
    sub_keys, sub_dict = key.split("/"), data
    for sub_key in sub_keys:
        sub_dict = sub_dict[sub_key]
    return sub_dict

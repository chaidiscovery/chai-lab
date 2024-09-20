# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

from typing import TypeVar

T = TypeVar("T")


def default(x: T | None, y: T) -> T:
    return x if x is not None else y

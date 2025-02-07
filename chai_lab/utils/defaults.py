# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from typing import TypeVar

T = TypeVar("T")


def default(x: T | None, y: T) -> T:
    return x if x is not None else y

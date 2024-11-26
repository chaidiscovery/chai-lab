# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import typing

from beartype import beartype
from jaxtyping import (
    Bool,
    Float,
    Float32,
    Int,
    Int32,
    Num,
    Shaped,
    TypeCheckError,
    UInt8,
    jaxtyped,
)

# Modules are only loaded and executed the first time they are imported, so the value of
# should_typecheck will constant over the lifetime of the program.
should_typecheck = True


Func = typing.TypeVar("Func")


def typecheck(cls_or_func: Func) -> Func:
    if should_typecheck:
        return jaxtyped(typechecker=beartype)(cls_or_func)  # type: ignore
    else:
        return cls_or_func


__all__ = [
    "typecheck",
    "TypeCheckError",
    # re-export jaxtyping types
    "Bool",
    "Float",
    "Int",
    "Int32",
    "Float32",
    "Num",
    "Shaped",
    "UInt8",
]

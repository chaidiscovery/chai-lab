from typing import TypeVar

T = TypeVar("T")


def default(x: T | None, y: T) -> T:
    return x if x is not None else y

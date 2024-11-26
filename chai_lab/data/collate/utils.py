# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)

# static graph is exported for different n_tokens,
#  we pad to the closest one
AVAILABLE_MODEL_SIZES = [256, 384, 512, 768, 1024, 1536, 2048]


@dataclass(frozen=True)
class PadSizes:
    n_tokens: int
    n_atoms: int


def pad_size(max_in_batch: int, allowed_sizes: list[int]) -> int:
    """pads to the smallest allowed size"""
    max_allowed_size = allowed_sizes[-1]
    if max_in_batch > max_allowed_size:
        raise ValueError(f"{max_in_batch=} > {max_allowed_size=}")
    return min(n for n in allowed_sizes if n >= max_in_batch)


def get_pad_sizes(contexts: list[AllAtomStructureContext]) -> PadSizes:
    max_n_tokens = max(context.num_tokens for context in contexts)
    n_tokens = pad_size(max_n_tokens, AVAILABLE_MODEL_SIZES)

    max_n_atoms = max(context.num_atoms for context in contexts)
    n_atoms = 23 * n_tokens
    assert max_n_atoms <= n_atoms
    assert n_atoms % 32 == 0

    return PadSizes(n_tokens=n_tokens, n_atoms=n_atoms)

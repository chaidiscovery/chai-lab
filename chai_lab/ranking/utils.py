# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from einops import rearrange
from torch import Tensor

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.tensor_utils import cdist
from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
def get_chain_masks_and_asyms(
    asym_id: Int[Tensor, "... n"],
    mask: Bool[Tensor, "... n"],
) -> tuple[Bool[Tensor, "... c n"], Int[Tensor, "c"]]:
    """
    Returns a mask for each chain and the unique asym ids
    """
    sorted_unique_asyms = torch.unique(asym_id[mask])
    # shape: (..., max_num_chains, n)
    chain_masks = rearrange(asym_id, "... n -> ... 1 n") == rearrange(
        sorted_unique_asyms, "nc -> nc 1"
    )  # shape: (..., n, max_num_chains)
    return chain_masks & rearrange(mask, "... n -> ... 1 n"), sorted_unique_asyms


@typecheck
def get_interface_mask(
    coords: Float[Tensor, "... n 3"],
    asym_id: Int[Tensor, "... n"],
    mask: Bool[Tensor, "... n"],
    interface_threshold: float,
) -> Bool[Tensor, "... n n"]:
    valid_mask = rearrange(asym_id, "... n -> ... n 1") != rearrange(
        asym_id, "... n -> ... 1 n"
    )
    valid_mask &= rearrange(mask, "... n -> ... n 1") & rearrange(
        mask, "... n -> ... 1 n"
    )
    dists = torch.masked_fill(cdist(coords), ~valid_mask, torch.inf)
    min_dists, _ = torch.min(dists, dim=-1)
    return min_dists < interface_threshold


@typecheck
def expectation(
    logits: Float[Tensor, "... bins"],
    weights: Float[Tensor, "... bins"],
) -> Float[Tensor, "..."]:  # last dim will be dropped
    logits = torch.softmax(logits, dim=-1)
    return (logits * weights).sum(dim=-1)


@typecheck
def num_atoms_per_chain(
    atom_mask: Bool[Tensor, "... a"],
    asym_id: Int[Tensor, "... a"],
) -> Int[Tensor, "... c"]:
    masks, _ = get_chain_masks_and_asyms(asym_id, atom_mask)
    return masks.sum(dim=-1)


@typecheck
def chain_is_polymer(
    asym_id: Int[Tensor, "... n"],
    mask: Bool[Tensor, "... n"],
    entity_type: Int[Tensor, "... n"],
) -> Bool[Tensor, "... c"]:
    chain_masks, _ = get_chain_masks_and_asyms(asym_id, mask)
    polymer_types = torch.tensor(
        [
            EntityType.PROTEIN.value,
            EntityType.RNA.value,
            EntityType.DNA.value,
            EntityType.POLYMER_HYBRID.value,
        ],
        device=entity_type.device,
    )
    is_polymer = torch.any(entity_type.unsqueeze(-1) == polymer_types, dim=-1)
    chain_is_polymer = []
    for polymer_mask in chain_masks.unbind(dim=-2):
        chain_is_polymer.append(torch.any(is_polymer & polymer_mask, dim=-1))
    return torch.stack(chain_is_polymer, dim=-1)

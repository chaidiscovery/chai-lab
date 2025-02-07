# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import torch
from einops import rearrange, reduce, repeat
from torch import Tensor

import chai_lab.ranking.utils as rank_utils
from chai_lab.utils.tensor_utils import cdist, und_self
from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
@dataclass
class ClashScores:
    """
    total_clashes: total number of clashes in the complex
    total_inter_chain_clashes: total number of inter-chain clashes in the complex,
        i.e. inter-chain clashes summed over all chain pairs
    per_chain_intra_clashes: number of intra-chain clashes for each chain in the complex
    per_chain_pair_clashes: number of inter-chain clashes for each chain pair in the complex
    """

    total_clashes: Int[Tensor, "..."]
    total_inter_chain_clashes: Int[Tensor, "..."]
    chain_chain_clashes: Int[Tensor, "... n_chains n_chains"]
    has_inter_chain_clashes: Bool[Tensor, "..."]


@typecheck
def _compute_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    clash_threshold: float = 1.1,
) -> Bool[Tensor, "... a a"]:
    pairwise_dists = cdist(atom_coords)
    valid_mask = und_self(atom_mask, "... i, ... j -> ... i j")
    valid_mask = valid_mask & ~torch.eye(
        atom_coords.shape[-2], device=atom_coords.device, dtype=torch.bool
    )
    return valid_mask & (pairwise_dists < clash_threshold)


@typecheck
def has_inter_chain_clashes(
    atom_mask: Bool[Tensor, "... a"],
    atom_asym_id: Int[Tensor, "... a"],
    atom_entity_type: Int[Tensor, "... a"],
    per_chain_pair_clashes: Int[Tensor, "... n_chains n_chains"],
    max_clashes: int = 100,
    max_clash_ratio: float = 0.5,
) -> Bool[Tensor, "..."]:
    """
    Determine if the complex has inter-chain clashes.
    Criteria:
        (1) If a chain pair has more than `max_clashes` clashes, then consider it a clash
        (2) If a chain pair has less than `max_clashes` clashes, but the total number of
            clashes is more than `max_clash_ratio` of the smaller chain's total atoms,
            then also consider it a clash
        (3) The chain pairs must be both be polymers

    """
    has_clashes = per_chain_pair_clashes >= max_clashes

    atoms_per_chain = rank_utils.num_atoms_per_chain(
        atom_mask=atom_mask,
        asym_id=atom_asym_id,
    )

    # if a chain pair has less than max_clashes clashes, butmore than
    # max_clash_ratio of the smaller chain's total atoms, then also
    # consider it a clash
    has_clashes |= (
        per_chain_pair_clashes
        / rearrange(atoms_per_chain, "... c -> ... c 1").clamp(min=1)
    ).ge(max_clash_ratio)

    has_clashes |= (
        per_chain_pair_clashes / rearrange(atoms_per_chain, "b c -> b 1 c").clamp(min=1)
    ).ge(max_clash_ratio)

    # only consider clashes between pairs of polymer chains
    polymer_chains = rank_utils.chain_is_polymer(
        asym_id=atom_asym_id,
        mask=atom_mask,
        entity_type=atom_entity_type,
    )
    is_polymer_pair = und_self(polymer_chains, "... c1, ... c2 -> ... c1 c2")

    # reduce over all chain pairs
    return reduce(has_clashes & is_polymer_pair, "... c1 c2 -> ...", torch.any)


@typecheck
def get_scores(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    atom_asym_id: Int[Tensor, "... a"],
    atom_entity_type: Int[Tensor, "... a"],
    clash_threshold: float = 1.1,
    max_clashes: int = 100,
    max_clash_ratio: float = 0.5,
) -> ClashScores:
    # shift asym_id from 1-based to 0-based
    assert atom_asym_id.dtype in (torch.int32, torch.int64)
    atom_asym_id = (atom_asym_id - 1).to(torch.int64)
    assert torch.amin(atom_asym_id) >= 0

    # dimensions
    n_chains = atom_asym_id.amax().add(1).item()
    assert isinstance(n_chains, int)
    *b, a = atom_mask.shape

    clashes_a_a = _compute_clashes(atom_coords, atom_mask, clash_threshold)
    clashes_a_a = clashes_a_a.to(torch.int32)  # b a a

    clashes_a_chain = clashes_a_a.new_zeros(*b, a, n_chains)
    clashes_a_chain.scatter_add_(
        dim=-1,
        src=clashes_a_a,
        index=repeat(atom_asym_id, f"b a -> b {a} a"),
    )

    clashes_chain_chain = clashes_a_a.new_zeros(*b, n_chains, n_chains)
    clashes_chain_chain.scatter_add_(
        dim=-2,
        src=clashes_a_chain,
        index=repeat(atom_asym_id, f"b a -> b a {n_chains}"),
    )
    # i, j enumerate chains
    total_clashes = reduce(clashes_chain_chain, "... i j -> ...", "sum") // 2

    # NB: self-interaction of chain contains doubled self-interaction,
    #  we compensate for this.
    clashes_chain_chain = clashes_chain_chain // (
        1 + torch.diag(clashes_a_a.new_ones(n_chains))
    )
    # in case anyone needs
    # per_chain_intra_clashes = torch.einsum("... i i -> ... i", clashes_chain_chain)
    # delete self-interaction for simplicity
    non_diag = 1 - torch.diag(clashes_a_a.new_ones(n_chains))
    inter_chain_chain = non_diag * clashes_chain_chain

    inter_chain_clashes = (
        reduce(inter_chain_chain, "... i j -> ... ", "sum") // 2
    )  # div by 2 to compensate for symmetricity of matrix

    return ClashScores(
        total_clashes=total_clashes,
        total_inter_chain_clashes=inter_chain_clashes,
        chain_chain_clashes=clashes_chain_chain,
        has_inter_chain_clashes=has_inter_chain_clashes(
            atom_mask=atom_mask,
            atom_asym_id=atom_asym_id,
            atom_entity_type=atom_entity_type,
            per_chain_pair_clashes=inter_chain_chain,
            max_clashes=max_clashes,
            max_clash_ratio=max_clash_ratio,
        ),
    )

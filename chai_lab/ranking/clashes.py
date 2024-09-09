from dataclasses import dataclass
from itertools import combinations

import torch
from einops import rearrange, reduce, repeat
from torch import Tensor

import chai_lab.ranking.utils as rutils
from chai_lab.ranking.utils import (
    get_chain_masks_and_asyms,
)
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

    total_clashes: Float[Tensor, "..."]
    total_inter_chain_clashes: Float[Tensor, "..."]
    per_chain_intra_clashes: Float[Tensor, "... n_chains"]
    per_chain_pair_clashes: Float[Tensor, "... n_chains n_chains"]
    has_clashes: Bool[Tensor, "..."]


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
def maybe_compute_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    clash_matrix: Bool[Tensor, "... a a"] | None = None,
    clash_threshold: float = 1.1,
) -> Bool[Tensor, "... a a"]:
    if clash_matrix is None:
        return _compute_clashes(atom_coords, atom_mask, clash_threshold)
    else:
        return clash_matrix


@typecheck
def total_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    clash_matrix: Bool[Tensor, "... a a"] | None = None,
    clash_threshold: float = 1.1,
) -> Float[Tensor, "..."]:
    """
    Computes the total number of clashes in the complex
    """
    clash_matrix = maybe_compute_clashes(
        atom_coords, atom_mask, clash_matrix, clash_threshold
    )
    # clash matrix is symmetric
    return reduce(clash_matrix, "... a1 a2 -> ...", "sum") / 2


@typecheck
def total_inter_chain_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    asym_id: Int[Tensor, "... a"],
    clash_matrix: Bool[Tensor, "... a a"] | None = None,
    clash_threshold: float = 1.1,
) -> Float[Tensor, "..."]:
    """Compute total number of inter-chain clashes in the complex"""
    clash_matrix = maybe_compute_clashes(
        atom_coords, atom_mask, clash_matrix, clash_threshold
    ).clone()  # don't overwrite an input
    # clash matrix is symmetric
    clash_matrix &= rearrange(asym_id, "... a -> ... a 1") != rearrange(
        asym_id, "... a -> ... 1 a"
    )
    # account for double counting
    return reduce(clash_matrix, "... a1 a2 -> ...", "sum") / 2


@typecheck
def per_chain_intra_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    asym_id: Int[Tensor, "... a"],
    clash_matrix: Bool[Tensor, "... a a"] | None = None,
    clash_threshold: float = 1.1,
) -> tuple[Float[Tensor, "... n_chains"], Int[Tensor, "n_chains"]]:
    clash_matrix = maybe_compute_clashes(
        atom_coords, atom_mask, clash_matrix, clash_threshold
    ).clone()  # don't overwrite an input
    # clash matrix is symmetric
    clash_matrix &= rearrange(asym_id, "... a -> ... a 1") == rearrange(
        asym_id, "... a -> ... 1 a"
    )
    per_atom_clashes = reduce(clash_matrix, "... a -> ...", "sum") / 2
    # add dimension for chains
    per_atom_clashes = rearrange(per_atom_clashes, "... a -> ... 1 a")
    chain_masks, asyms = get_chain_masks_and_asyms(asym_id, atom_mask)
    return reduce(per_atom_clashes * chain_masks, "... c a -> ... c", "sum"), asyms


@typecheck
def per_chain_pair_clashes(
    atom_coords: Float[Tensor, "... a 3"],
    atom_mask: Bool[Tensor, "... a"],
    asym_id: Int[Tensor, "... a"],
    clash_matrix: Bool[Tensor, "... a a"] | None = None,
    clash_threshold: float = 1.1,
) -> tuple[Float[Tensor, "... n_chains n_chains"], Int[Tensor, "n_chains"]]:
    """
    Compute the number of inter-chain clashes for each chain in the complex
    """
    clash_matrix = maybe_compute_clashes(
        atom_coords, atom_mask, clash_matrix, clash_threshold
    ).clone()  # don't overwrite an input
    clash_matrix &= rearrange(asym_id, "... a -> ... a 1") != rearrange(
        asym_id, "... a -> ... 1 a"
    )
    chain_masks, asyms = get_chain_masks_and_asyms(asym_id, atom_mask)
    per_chain_clashes = torch.zeros(
        *chain_masks.shape[:-2],
        len(asyms),
        len(asyms),
        device=atom_coords.device,
        dtype=torch.float32,
    )
    # compute in loop to minimize peak memory
    for i, j in combinations(range(len(asyms)), 2):
        chain_pair_mask = torch.einsum(
            "...i,...j->...ij", chain_masks[..., i, :], chain_masks[..., j, :]
        )
        # chain_pair_mask is triangular, so don't need to account for double counting
        per_chain_clashes[..., i, j] = reduce(
            clash_matrix * chain_pair_mask, "... i j -> ...", "sum"
        )
    symm_clashes = per_chain_clashes + rearrange(
        per_chain_clashes, "... i j -> ... j i"
    )
    return symm_clashes, asyms


@typecheck
def has_clashes(
    atom_mask: Bool[Tensor, "... a"],
    atom_asym_id: Int[Tensor, "... a"],
    atom_entity_type: Int[Tensor, "... a"],
    per_chain_pair_clashes: Float[Tensor, "... n_chains n_chains"],
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

    atoms_per_chain = rutils.num_atoms_per_chain(
        atom_mask=atom_mask,
        asym_id=atom_asym_id,
    )

    # if a chain pair has less than max_clashes clashes, butmore than
    # max_clash_ratio of the smaller chain's total atoms, then also
    # consider it a clash
    c = atoms_per_chain.shape[-1]
    atoms_per_chain_row = repeat(atoms_per_chain, "... c -> ... (c k)", k=c)
    atoms_per_chain_col = repeat(atoms_per_chain, "... c -> ... (k c)", k=c)
    min_atoms_per_chain_pair, _ = torch.min(
        torch.stack([atoms_per_chain_row, atoms_per_chain_col], dim=-1), dim=-1
    )
    min_atoms_per_chain_pair = rearrange(
        min_atoms_per_chain_pair,
        "... (c_row c_col) -> ... c_row c_col",
        c_row=c,
    )
    has_clashes |= (
        per_chain_pair_clashes / torch.clamp(min_atoms_per_chain_pair, min=1)
    ) >= max_clash_ratio

    # only consider clashes between pairs of polymer chains
    polymer_chains = rutils.chain_is_polymer(
        asym_id=atom_asym_id,
        mask=atom_mask,
        entity_type=atom_entity_type,
    )
    is_polymer_pair = rearrange(polymer_chains, "... c -> ... c 1") & rearrange(
        polymer_chains, "... c -> ... 1 c"
    )
    # reduce over all chain pairs
    return torch.any(has_clashes & is_polymer_pair, dim=(-1, -2))


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
    clash_matrix = _compute_clashes(atom_coords, atom_mask, clash_threshold)
    _per_chain_pair_clashes = per_chain_pair_clashes(
        atom_coords, atom_mask, atom_asym_id, clash_matrix, clash_threshold
    )[0]
    return ClashScores(
        total_clashes=total_clashes(
            atom_coords=atom_coords,
            atom_mask=atom_mask,
            clash_matrix=clash_matrix,
            clash_threshold=clash_threshold,
        ),
        total_inter_chain_clashes=total_inter_chain_clashes(
            atom_coords=atom_coords,
            atom_mask=atom_mask,
            asym_id=atom_asym_id,
            clash_matrix=clash_matrix,
            clash_threshold=clash_threshold,
        ),
        per_chain_intra_clashes=per_chain_intra_clashes(
            atom_coords=atom_coords,
            atom_mask=atom_mask,
            asym_id=atom_asym_id,
            clash_matrix=clash_matrix,
            clash_threshold=clash_threshold,
        )[0],
        per_chain_pair_clashes=_per_chain_pair_clashes,
        has_clashes=has_clashes(
            atom_mask=atom_mask,
            atom_asym_id=atom_asym_id,
            atom_entity_type=atom_entity_type,
            per_chain_pair_clashes=_per_chain_pair_clashes,
            max_clashes=max_clashes,
            max_clash_ratio=max_clash_ratio,
        ),
    )

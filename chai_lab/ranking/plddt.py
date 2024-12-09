# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

from einops import repeat
from torch import Tensor

import chai_lab.ranking.utils as rank_utils
from chai_lab.utils.tensor_utils import masked_mean
from chai_lab.utils.typing import Bool, Float, Int, typecheck


@typecheck
@dataclass
class PLDDTScores:
    """
    complex_plddt: plddt score of the complex
    per_chain_plddt: plddt score for each chain in the complex
    per_atom_plddt: plddt score for each atom in the complex
    """

    complex_plddt: Float[Tensor, "..."]
    per_chain_plddt: Float[Tensor, "... c"]
    per_atom_plddt: Float[Tensor, "... a"]


@typecheck
def plddt(
    logits: Float[Tensor, "... a bins"],
    mask: Bool[Tensor, "... a"],
    bin_centers: Float[Tensor, "bins"],
    per_residue: bool = False,
) -> Float[Tensor, "..."] | Float[Tensor, "... a"]:
    expectations = rank_utils.expectation(logits, bin_centers)
    if per_residue:
        return expectations
    else:
        return masked_mean(mask, expectations, dim=-1)


@typecheck
def per_chain_plddt(
    logits: Float[Tensor, "... a bins"],
    atom_mask: Bool[Tensor, "... a"],
    asym_id: Int[Tensor, "... a"],
    bin_centers: Float[Tensor, "bins"],
) -> Float[Tensor, "... c"]:
    chain_masks, _ = rank_utils.get_chain_masks_and_asyms(asym_id, atom_mask)
    logits = repeat(logits, "... a b -> ... c a b", c=chain_masks.shape[-2])
    return plddt(logits, chain_masks, bin_centers, per_residue=False)


@typecheck
def get_scores(
    lddt_logits: Float[Tensor, "... a bins"],
    atom_mask: Bool[Tensor, "... a"],
    atom_asym_id: Int[Tensor, "... a"],
    bin_centers: Float[Tensor, "bins"],
) -> PLDDTScores:
    return PLDDTScores(
        complex_plddt=plddt(
            logits=lddt_logits,
            mask=atom_mask,
            bin_centers=bin_centers,
            per_residue=False,
        ),
        per_atom_plddt=plddt(
            logits=lddt_logits,
            mask=atom_mask,
            bin_centers=bin_centers,
            per_residue=True,
        ),
        per_chain_plddt=per_chain_plddt(
            logits=lddt_logits,
            atom_mask=atom_mask,
            asym_id=atom_asym_id,
            bin_centers=bin_centers,
        ),
    )

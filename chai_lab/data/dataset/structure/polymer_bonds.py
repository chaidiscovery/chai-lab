# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""Peptide-link inference for modified CCD residues."""

import json
import logging
from itertools import pairwise
from pathlib import Path
from typing import TypedDict, cast

import torch

from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import standard_residue_pdb_codes

logger = logging.getLogger(__name__)


class PolymerLinkage(TypedDict):
    """Atoms used to join a CCD component into a peptide."""

    n_atom: str | None
    c_atom: str | None
    n_leaving_atoms: list[str]
    c_leaving_atoms: list[str]


_POLYMER_LINKAGES = cast(
    dict[str, PolymerLinkage],
    json.loads(Path(__file__).with_name("polymer_linkages.json").read_text()),
)
_STANDARD_LINKAGE = PolymerLinkage(
    n_atom="N",
    c_atom="C",
    n_leaving_atoms=[],
    c_leaving_atoms=[],
)


def add_inferred_polymer_bonds(context: AllAtomStructureContext) -> None:
    """Add peptide bonds involving modified residues and mask leaving atoms."""
    residues_by_chain: dict[int, dict[int, str]] = {}
    for token_index in range(context.num_tokens):
        if context.token_entity_type[token_index].item() != EntityType.PROTEIN.value:
            continue
        asym_id = context.token_asym_id[token_index].item()
        residue_index = context.token_residue_index[token_index].item()
        residue_name = context.residue_names[token_index]
        residues = residues_by_chain.setdefault(asym_id, {})
        previous_name = residues.setdefault(residue_index, residue_name)
        if previous_name != residue_name:
            raise ValueError(
                f"Residue {residue_index + 1} has conflicting names: "
                f"{previous_name} and {residue_name}"
            )

    atom_lookup: dict[tuple[int, int, str], int] = {}
    for atom_index, (token_index, atom_name) in enumerate(
        zip(context.atom_token_index.tolist(), context.atom_ref_name, strict=True)
    ):
        if context.token_entity_type[token_index].item() != EntityType.PROTEIN.value:
            continue
        atom_key = (
            context.token_asym_id[token_index].item(),
            context.token_residue_index[token_index].item(),
            atom_name,
        )
        if atom_key in atom_lookup:
            raise ValueError(f"Duplicate atom name in residue: {atom_key}")
        atom_lookup[atom_key] = atom_index

    inferred: list[tuple[int, int]] = []
    leaving_atoms: set[int] = set()
    for asym_id, residues in residues_by_chain.items():
        residue_indices = sorted(residues)
        for left_index, right_index in pairwise(residue_indices):
            if right_index != left_index + 1:
                continue
            left_name = residues[left_index]
            right_name = residues[right_index]
            if (
                left_name in standard_residue_pdb_codes
                and right_name in standard_residue_pdb_codes
            ):
                continue

            left_linkage = (
                _STANDARD_LINKAGE
                if left_name in standard_residue_pdb_codes
                else _POLYMER_LINKAGES.get(left_name)
            )
            right_linkage = (
                _STANDARD_LINKAGE
                if right_name in standard_residue_pdb_codes
                else _POLYMER_LINKAGES.get(right_name)
            )
            if (
                left_linkage is None
                or right_linkage is None
                or left_linkage["c_atom"] is None
                or right_linkage["n_atom"] is None
            ):
                logger.warning(
                    "Cannot infer peptide bond between %s %d and %s %d",
                    left_name,
                    left_index + 1,
                    right_name,
                    right_index + 1,
                )
                continue

            left_atom = atom_lookup.get((asym_id, left_index, left_linkage["c_atom"]))
            right_atom = atom_lookup.get(
                (asym_id, right_index, right_linkage["n_atom"])
            )
            if left_atom is None or right_atom is None:
                logger.warning(
                    "Missing peptide-link atoms between %s %d and %s %d",
                    left_name,
                    left_index + 1,
                    right_name,
                    right_index + 1,
                )
                continue

            inferred.extend([(left_atom, right_atom), (right_atom, left_atom)])
            for atom_name in left_linkage["c_leaving_atoms"]:
                atom_index = atom_lookup.get((asym_id, left_index, atom_name))
                if atom_index is not None:
                    leaving_atoms.add(atom_index)
            for atom_name in right_linkage["n_leaving_atoms"]:
                atom_index = atom_lookup.get((asym_id, right_index, atom_name))
                if atom_index is not None:
                    leaving_atoms.add(atom_index)

    existing = list(
        zip(
            context.atom_covalent_bond_indices[0].tolist(),
            context.atom_covalent_bond_indices[1].tolist(),
            strict=True,
        )
    )
    bonds = list(dict.fromkeys([*existing, *inferred]))
    if bonds:
        device = context.atom_covalent_bond_indices[0].device
        context.atom_covalent_bond_indices = (
            torch.tensor([bond[0] for bond in bonds], dtype=torch.long, device=device),
            torch.tensor([bond[1] for bond in bonds], dtype=torch.long, device=device),
        )
    if leaving_atoms:
        context.atom_exists_mask[list(leaving_atoms)] = False

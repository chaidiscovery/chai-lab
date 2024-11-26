# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from functools import lru_cache

import torch
from torch import Tensor

import chai_lab.data.residue_constants as rc
from chai_lab.utils.typing import Bool, Int


def get_centre_atom_name(residue_name: str) -> str:
    if residue_name not in rc.standard_residue_pdb_codes:
        raise ValueError(f"Residue {residue_name} is not a standard residue")

    if residue_name in {
        "A",
        "G",
        "C",
        "U",
        "DA",
        "DG",
        "DC",
        "DT",
    }:
        return "C1'"
    else:
        assert len(residue_name) == 3, "residue expected"
        return "CA"


def get_reference_atom_name(residue_name: str) -> str:
    if residue_name not in rc.standard_residue_pdb_codes:
        raise ValueError(f"Residue {residue_name} is not a standard residue")

    if residue_name == "GLY":
        return "CA"
    elif residue_name in {"A", "G", "DA", "DG"}:
        return "C4"
    elif residue_name in {"C", "U", "DC", "DT"}:
        return "C2"
    else:
        return "CB"


def get_centre_atom_index(atom_names: list[str], residue_name: str) -> Int[Tensor, "1"]:
    # centre of the token is  Calpha or C1'
    name = get_centre_atom_name(residue_name)

    if name in atom_names:
        idx = atom_names.index(name)
    else:
        raise ValueError(
            f"Residue {residue_name} marked as standard, "
            f"but reference conformer misses centre atom {name}. "
            "Either the residue is not standard or reference conformer is wrong."
        )

    return torch.tensor([idx], dtype=torch.int)


def get_reference_atom_index(
    atom_names: list[str], residue_name: str
) -> Int[Tensor, "1"]:
    name = get_reference_atom_name(residue_name)
    if name in atom_names:
        idx = atom_names.index(name)
    else:
        raise ValueError(
            f"Residue {residue_name} marked as standard, "
            f"but reference conformer misses reference atom {name}. "
            "Either the residue is not standard or reference conformer is wrong."
        )

    return torch.tensor([idx], dtype=torch.int)


def get_backbone_frame_atom_names(residue_name: str) -> tuple[str, str, str]:
    """Return names of the 3 atoms used in canonical token frame."""
    if residue_name in {
        "A",
        "G",
        "C",
        "U",
        "DA",
        "DG",
        "DC",
        "DT",
    }:
        return "C1'", "C3'", "C4'"
    if residue_name in rc.residue_atoms:
        return "N", "CA", "C"
    return "", "", ""


def backbone_atoms_all_present(
    atom_names: list[str], residue_name: str
) -> Bool[Tensor, "1"]:
    """Check if all *protein* backbone atoms are present in the list of atom names."""
    backbone_frame_atoms = get_backbone_frame_atom_names(residue_name)
    if all(a == "" for a in backbone_frame_atoms):
        # Not a nucleic acid or a protein residue
        all_present = False
    else:
        all_present = all(name in atom_names for name in backbone_frame_atoms)
    return torch.tensor([all_present], dtype=torch.bool)


def backbone_atoms_indices(
    atom_names: list[str], residue_name: str
) -> Int[Tensor, "3"]:
    """Return indices of backbone atoms N, Ca, C in the list of atom names."""
    backbone_frame_atom_names = get_backbone_frame_atom_names(residue_name)

    if backbone_atoms_all_present(atom_names, residue_name):
        indices = [atom_names.index(name) for name in backbone_frame_atom_names]
    else:
        indices = [0, 0, 0]

    return torch.tensor(indices, dtype=torch.int)


@lru_cache(maxsize=1)
def atom_37_atom_indices() -> dict[tuple[str, str | None], int]:
    num_protein_atoms = 37
    protein_res_atom_to_index: dict[tuple[str, str | None], int] = {
        (residue_name, atom_name): atom_index
        for residue_name in rc.residue_atoms.keys()
        for atom_name, atom_index in rc.atom_order.items()
    }
    assert max(protein_res_atom_to_index.values()) == num_protein_atoms - 1

    num_rna_atoms = 36
    # note: convert RNA residues to R{} to match residue names from residue_constants.py
    rna_res_atom_to_index = {
        (residue_name, atom_name): atom_index
        for residue_name in {"A", "C", "G", "U"}
        for atom_index, atom_name in enumerate(
            rc.nucleic_acid_atoms[f"R{residue_name}"]
        )
    }
    assert max(rna_res_atom_to_index.values()) == num_rna_atoms - 1

    num_dna_atoms = 36
    dna_res_atom_to_index = {
        (residue_name, atom_name): atom_index
        for residue_name in {"DA", "DC", "DG", "DT"}
        for atom_index, atom_name in enumerate(rc.nucleic_acid_atoms[residue_name])
    }
    assert max(dna_res_atom_to_index.values()) == num_dna_atoms - 1

    return {
        **protein_res_atom_to_index,
        **rna_res_atom_to_index,
        **dna_res_atom_to_index,
    }

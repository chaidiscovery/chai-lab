# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import gemmi
import numpy as np
import torch
from torch import Tensor

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import residue_types_with_nucleotides_order
from chai_lab.model.utils import center_random_augmentation
from chai_lab.utils.typing import Bool, Float, Int


@dataclass
class ConformerData:
    position: Float[Tensor, "n 3"]
    element: Int[Tensor, "n"]
    charge: Int[Tensor, "n"]
    atom_names: list[str]
    bonds: list[tuple[int, int]]
    symmetries: Int[Tensor, "n n_symm"]

    @property
    def num_atoms(self) -> int:
        num_atoms, _ = self.position.shape
        assert num_atoms == len(self.atom_names)
        return num_atoms

    def gather_atom_positions(
        self, query_atom_names: list[str]
    ) -> tuple[Float[Tensor, "n 3"], Bool[Tensor, "n"]]:
        if self.num_atoms == 0:
            gathered_positions = torch.zeros(len(query_atom_names), 3)
            mask = torch.zeros(len(query_atom_names), dtype=torch.bool)
            return gathered_positions, mask

        atom_indices = {name: i for i, name in enumerate(self.atom_names)}
        indices = torch.tensor(
            [atom_indices.get(name, -1) for name in query_atom_names],
            dtype=torch.int,
        )
        mask = indices != -1
        gathered_positions = self.position[indices] * mask.unsqueeze(-1)

        return gathered_positions, mask

    def center_random_augment(
        self,
    ) -> "ConformerData":
        if self.num_atoms == 0:
            return self

        atom_mask = torch.ones_like(self.element, dtype=torch.bool)
        centered_coords = center_random_augmentation(
            self.position.unsqueeze(0), atom_mask.unsqueeze(0)
        )[0]
        return ConformerData(
            centered_coords,
            self.element,
            self.charge,
            self.atom_names,
            self.bonds,
            self.symmetries,
        )


@dataclass
class Residue:
    name: str
    label_seq: int | None
    restype: int
    residue_index: int
    is_missing: bool
    b_factor_or_plddt: float
    conformer_data: ConformerData | None
    smiles: str | None = None
    is_covalent_bonded: bool = False


def get_restype(
    residue_info: gemmi.ResidueInfo,
    entity_type: EntityType,
) -> int:
    """
    Encodes residues into alphabet of size 32:
    20 standards AAs + X
    4 RNA bases + RX
    4 DNA bases + DX
    GAP
    note: ligand residues as encoded as X
    """

    if residue_info.is_amino_acid():
        restype = residue_info.fasta_code()  # encodes non-standard as X
        unknown_value = residue_types_with_nucleotides_order["X"]
    elif residue_info.is_nucleic_acid() and entity_type == EntityType.RNA:
        restype = "R{}".format(residue_info.one_letter_code)
        unknown_value = residue_types_with_nucleotides_order["RX"]
    elif residue_info.is_nucleic_acid() and entity_type == EntityType.DNA:
        restype = "D{}".format(residue_info.one_letter_code)
        unknown_value = residue_types_with_nucleotides_order["DX"]
    else:
        restype = "X"
        unknown_value = residue_types_with_nucleotides_order["X"]

    tokenized_restype = residue_types_with_nucleotides_order.get(restype, unknown_value)
    return tokenized_restype


def _missing_residue(label_seq: int, name: str) -> gemmi.Residue:
    """
    Creates a missing residue with given type (three letter code)
    and label_seq (1-indexed position in the reference protein sequence for given entity)
    """
    res = gemmi.Residue()
    res.label_seq = label_seq
    res.name = name
    return res


def _get_residue(
    residue: gemmi.Residue,
    entity_type: EntityType,
    is_covalent_bonded: bool = False,
) -> Residue:
    residue_info = gemmi.find_tabulated_residue(residue.name)
    restype = get_restype(residue_info, entity_type)

    # note: residue.label_seq is the residue index
    # in the full sequence of each chain that was crystallized
    # it may not be contiguous if the crystallized structure has missing residues
    # label_seq is 1-indexed, we make a 0-indexed residue_index
    residue_index = 0 if residue.label_seq is None else residue.label_seq - 1

    is_missing = len([atom for atom in residue.first_conformer()]) == 0

    # The ground truth conformer is the data that we obtain from the PDB structure.
    conformer_data = _get_gemmi_conformer_data(residue)

    b_factor_or_plddt = (
        np.mean([atom.b_iso for atom in residue.first_conformer()]).item()
        if not is_missing
        else 0.0
    )

    return Residue(
        name=residue.name,
        label_seq=residue.label_seq,
        restype=restype,
        residue_index=residue_index,
        is_missing=is_missing,
        b_factor_or_plddt=b_factor_or_plddt,
        conformer_data=conformer_data,
        is_covalent_bonded=is_covalent_bonded,
    )


def _select_among_altlocs(residue: gemmi.Residue, atom_name: str):
    # get one with largest occupancy
    group = residue[atom_name]
    # Max over a list of (occupancy, -index) -> extract the index
    # The -i means we tie-break by the first occurrence/lowest index
    atom_index = max([(a.occ, -i) for i, a in enumerate(group)])[1]
    return group[atom_index]


def get_heavy_atoms(residue: gemmi.Residue) -> dict[str, gemmi.Atom]:
    heavy_atoms = [
        atom for atom in residue.first_conformer() if atom.element.name != "H"
    ]

    atom_names = {atom.name for atom in heavy_atoms}

    return {name: _select_among_altlocs(residue, name) for name in atom_names}


def _get_gemmi_conformer_data(
    residue: gemmi.Residue,
) -> ConformerData:
    # returns a conformer data where atoms are sorted alphabetically by atom name
    # for reproducibility
    atoms = get_heavy_atoms(residue)
    ordered_atom_names = sorted(atoms.keys())

    if len(atoms) == 0:
        return ConformerData(
            position=torch.zeros(0, 3),
            charge=torch.zeros(0, dtype=torch.int),
            element=torch.zeros(0, dtype=torch.int),
            atom_names=[],
            bonds=[],
            symmetries=torch.zeros(0, 0),
        )

    position = torch.stack(
        [torch.tensor(atoms[k].pos.tolist()) for k in ordered_atom_names],
        dim=0,
    )
    charge = torch.tensor(
        [atoms[k].charge for k in ordered_atom_names],
        dtype=torch.int,
    )

    element = torch.tensor(
        [atoms[k].element.atomic_number for k in ordered_atom_names],
        dtype=torch.int,
    )

    return ConformerData(
        position=position,
        charge=charge,
        element=element,
        atom_names=ordered_atom_names,
        bonds=[],
        symmetries=torch.arange(len(atoms)).unsqueeze(-1),
    )


def get_residues(
    subchain: gemmi.ResidueSpan,
    full_sequence: list[str],
    entity_type: EntityType,
) -> list[Residue]:
    # Create a mapping from the label seq to the resolved entities.
    resolved_residues = {
        residue.label_seq: residue for residue in subchain.first_conformer()
    }

    # Make sure we add all residues from entity full sequence:
    all_residues = [
        resolved_residues.get(
            label_seq, _missing_residue(label_seq=label_seq, name=res_name)
        )
        for label_seq, res_name in enumerate(full_sequence, start=1)
    ]

    return [
        _get_residue(residue=residue, entity_type=entity_type)
        for residue in all_residues
    ]

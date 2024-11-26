# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import dataclass

import gemmi
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

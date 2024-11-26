# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
import string
from collections import defaultdict
from dataclasses import asdict, dataclass
from functools import cached_property
from pathlib import Path

import gemmi
from torch import Tensor

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.tensor_utils import tensorcode_to_string
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


_CHAIN_VOCAB = string.ascii_uppercase + string.ascii_lowercase


def get_pdb_chain_name(asym_id: int) -> str:
    vocab_index = asym_id - 1  # 1 -> A
    if vocab_index >= len(_CHAIN_VOCAB):
        # two letter codes
        chr1, chr2 = (
            (vocab_index // len(_CHAIN_VOCAB)) - 1,
            vocab_index % len(_CHAIN_VOCAB),
        )
        return _CHAIN_VOCAB[chr1] + _CHAIN_VOCAB[chr2]
    return _CHAIN_VOCAB[vocab_index]


@dataclass(frozen=True)
class PDBAtom:
    record_type: str
    atom_index: int
    atom_name: str
    alt_loc: str
    res_name_3: str
    chain_tag: str
    asym_id: int
    residue_index: int
    insertion_code: str
    pos: list[float]
    occupancy: float
    b_factor: float
    element: str
    charge: str

    def __str__(
        self,
    ):
        # currently this works only for single-char chain tags
        atom_line = (
            f"{self.record_type:<6}{self.atom_index:>5} {self.atom_name:<4}{self.alt_loc:>1}"
            f"{self.res_name_3:>3} {self.chain_tag:>1}"
            f"{self.residue_index:>4}{self.insertion_code:>1}   "
            f"{self.pos[0]:>8.3f}{self.pos[1]:>8.3f}{self.pos[2]:>8.3f}"
            f"{self.occupancy:>6.2f}{self.b_factor:>6.2f}          "
            f"{self.element:>2}{self.charge:>2}"
        )
        return atom_line

    def rename(self, atom_name: str) -> "PDBAtom":
        return PDBAtom(
            self.record_type,
            self.atom_index,
            atom_name,
            self.alt_loc,
            self.res_name_3,
            self.chain_tag,
            self.asym_id,
            self.residue_index,
            self.insertion_code,
            self.pos,
            self.occupancy,
            self.b_factor,
            self.element,
            self.charge,
        )


def write_pdb(chain_atoms: list[list[PDBAtom]], out_path: str):
    with open(out_path, "w") as f:
        for chain in chain_atoms:
            for atom in chain:
                f.write(str(atom) + "\n")
            f.write("TER\n")
        f.write("END\n")


@typecheck
@dataclass
class PDBContext:
    """Complex (multiple entities) represented as a collection of tensors"""

    token_residue_index: Int[Tensor, "n_tokens"]
    token_asym_id: Int[Tensor, "n_tokens"]
    token_entity_type: Int[Tensor, "n_tokens"]
    token_entity_id: Int[Tensor, "n_tokens"]
    token_residue_names: UInt8[Tensor, "n_tokens 8"]
    token_centre_atom_index: Int[Tensor, "n_tokens"]
    atom_token_index: Int[Tensor, "n_atoms"]
    atom_ref_element: Int[Tensor, "n_atoms"]
    atom_ref_mask: Bool[Tensor, "n_atoms"]
    atom_coords: Float[Tensor, "n_atoms 3"]
    atom_exists_mask: Bool[Tensor, "n_atoms"]
    token_exists_mask: Bool[Tensor, "n_tokens"]
    atom_ref_name_chars: Int[Tensor, "n_atoms 4"]
    atom_within_token_index: Int[Tensor, "n_atoms"]
    atom_bfactor_or_plddt: Float[Tensor, "n_atoms"] | None = None

    @cached_property
    def token_res_names_to_string(self) -> list[str]:
        return [tensorcode_to_string(x) for x in self.token_residue_names.cpu()]

    @property
    def is_ligand(self) -> bool:
        return self.is_entity(EntityType.LIGAND)

    def is_entity(self, ety: EntityType) -> bool:
        return self.token_entity_type[0].item() == ety.value

    def get_chain_entity_type(self, asym_id: int) -> int:
        mask = self.token_asym_id == asym_id
        assert mask.sum() > 0
        e_type = self.token_entity_type[mask][0].item()
        assert isinstance(e_type, int)
        return e_type

    def get_pdb_atoms(self):
        # warning: calling this on cuda tensors is extremely slow
        atom_asym_id = self.token_asym_id[self.atom_token_index]
        # atom level attributes
        atom_residue_index = (
            self.token_residue_index[self.atom_token_index] + 1
        )  # residues are 1-indexed
        atom_names = _tensor_to_atom_names(self.atom_ref_name_chars)
        atom_res_names = self.token_residue_names[self.atom_token_index]
        atom_res_names_strs = [
            tensorcode_to_string(x)[:3].ljust(3) for x in atom_res_names
        ]
        atom_element_names = [
            _atomic_num_to_element(int(x.item())) for x in self.atom_ref_element
        ]

        pdb_atoms = []
        num_atoms = self.atom_coords.shape[0]
        for atom_index in range(num_atoms):
            if not self.atom_exists_mask[atom_index].item():
                # skip missing atoms
                continue

            atom = PDBAtom(
                record_type="ATOM" if not self.is_ligand else "HETATM",
                atom_index=atom_index,
                atom_name=atom_names[atom_index],
                alt_loc="",
                res_name_3=atom_res_names_strs[atom_index],
                chain_tag=get_pdb_chain_name(int(atom_asym_id[atom_index].item())),
                asym_id=int(atom_asym_id[atom_index].item()),
                residue_index=int(atom_residue_index[atom_index].item()),
                insertion_code="",
                pos=self.atom_coords[atom_index].tolist(),
                occupancy=1.00,
                b_factor=(
                    1.00
                    if self.atom_bfactor_or_plddt is None
                    else self.atom_bfactor_or_plddt[atom_index].item()
                ),
                element=atom_element_names[atom_index],
                charge="",
            )
            pdb_atoms.append(atom)
        return pdb_atoms


def _atomic_num_to_element(atomic_num: int) -> str:
    return gemmi.Element(atomic_num).name


def entity_to_pdb_atoms(entity: PDBContext) -> list[list[PDBAtom]]:
    """Writes a single tokenized entity to PDB file"""
    pdb_atoms = entity.get_pdb_atoms()
    chains = defaultdict(list)
    for atom in pdb_atoms:
        chains[atom.asym_id].append(atom)

    for asym_id in chains:
        # deduplicate atom names within all ligand chains
        if entity.get_chain_entity_type(asym_id) == EntityType.LIGAND.value:
            chains[asym_id] = rename_ligand_atoms(chains[asym_id])

    return list(chains.values())


def rename_ligand_atoms(atoms: list[PDBAtom]) -> list[PDBAtom]:
    # transform ligand atom names from "C" to "C1", "C2"...
    atom_type_counter: dict[str, int] = {}
    renumbered_atoms = []
    for atom in atoms:
        idx = atom_type_counter.get(atom.element, 1)
        atom_type_counter[atom.element] = idx + 1
        base_name = atom.atom_name
        renumbered_atoms.append(atom.rename(f"{base_name}_{idx}"))
    return renumbered_atoms


def entities_to_pdb_file(entities: list[PDBContext], path: str):
    pdb_atoms: list[list[PDBAtom]] = []
    for entity in entities:
        pdb_atoms = pdb_atoms + entity_to_pdb_atoms(entity)
    write_pdb(pdb_atoms, path)


def pdb_context_from_batch(
    d: dict, coords: Tensor, plddt: Tensor | None = None
) -> PDBContext:
    context = PDBContext(
        token_residue_index=d["token_residue_index"][0],
        token_asym_id=d["token_asym_id"][0],
        token_entity_type=d["token_entity_type"][0],
        token_entity_id=d["token_entity_id"][0],
        token_residue_names=d["token_residue_name"][0],
        token_centre_atom_index=d["token_centre_atom_index"][0],
        atom_token_index=d["atom_token_index"][0],
        atom_ref_element=d["atom_ref_element"][0],
        atom_ref_mask=d["atom_ref_mask"][0],
        atom_coords=coords[0].cpu(),
        atom_exists_mask=d["atom_exists_mask"][0],
        token_exists_mask=d["token_exists_mask"][0],
        atom_ref_name_chars=d["atom_ref_name_chars"][0],
        atom_bfactor_or_plddt=plddt[0].cpu() if plddt is not None else None,
        atom_within_token_index=d["atom_within_token_index"][0],
    )
    for k, v in asdict(context).items():
        assert v.device.type == "cpu", ("not on cpu:", k, v.device)
    return context


def write_pdbs_from_outputs(
    coords: Float[Tensor, "1 n_atoms 3"],
    output_batch: dict,
    write_path: Path,
    bfactors: Float[Tensor, "1 n_atoms"] | None = None,
):
    # save outputs
    context = pdb_context_from_batch(output_batch, coords, plddt=bfactors)
    write_path.parent.mkdir(parents=True, exist_ok=True)
    entities_to_pdb_file(
        [context],
        str(write_path),
    )
    logger.info(f"saved pdb file to {write_path}")


@typecheck
def _tensor_to_atom_names(
    tensor: Int[Tensor, "*dims 4"] | UInt8[Tensor, "*dims 4"],
) -> list[str]:
    return [
        "".join([chr(ord_val + 32) for ord_val in ords_atom]).rstrip()
        for ords_atom in tensor
    ]

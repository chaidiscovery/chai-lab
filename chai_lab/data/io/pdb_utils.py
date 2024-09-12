import logging
import string
from collections import defaultdict
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import gemmi
from torch import Tensor

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.tensor_utils import tensorcode_to_string
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class PDBAtom:
    record_type: str
    atom_index: int
    atom_name: str
    alt_loc: str
    res_name_3: str
    chain_tag: str
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
    """Data needed to produce Posebuster input file types"""

    token_residue_index: Int[Tensor, "n_tokens"]
    token_asym_id: Int[Tensor, "n_tokens"]
    token_entity_type: Int[Tensor, "n_tokens"]
    token_residue_names: UInt8[Tensor, "n_tokens 8"]
    atom_token_index: Int[Tensor, "n_atoms"]
    atom_ref_element: Int[Tensor, "n_atoms"]
    atom_ref_mask: Bool[Tensor, "n_atoms"]
    atom_coords: Float[Tensor, "n_atoms 3"]
    atom_exists_mask: Bool[Tensor, "n_atoms"]
    atom_ref_name_chars: Int[Tensor, "n_atoms 4"]
    atom_bfactor_or_plddt: Float[Tensor, "n_atoms"] | None = None

    @cached_property
    def token_res_names_to_string(self) -> list[str]:
        return [tensorcode_to_string(x) for x in self.token_residue_names.cpu()]

    @property
    def num_atoms(self) -> int:
        return self.atom_coords.shape[0]

    @property
    def is_protein(self) -> bool:
        return self.is_entity(EntityType.PROTEIN)

    @property
    def is_ligand(self) -> bool:
        return self.is_entity(EntityType.LIGAND)

    @property
    def first_residue_name(self) -> str:
        return self.token_res_names_to_string[0].strip()

    def is_entity(self, ety: EntityType) -> bool:
        return self.token_entity_type[0].item() == ety.value

    def get_pdb_atoms(self):
        # warning: calling this on cuda tensors is extremely slow
        atom_asym_id = self.token_asym_id[self.atom_token_index]
        # atom level attributes
        atom_residue_index = (
            self.token_residue_index[self.atom_token_index] + 1
        )  # residues are 1-indexed
        atom_names = _tensor_to_atom_names(self.atom_ref_name_chars.unsqueeze(0))
        atom_res_names = self.token_residue_names[self.atom_token_index]
        atom_res_names_strs = [
            tensorcode_to_string(x)[:3].ljust(3) for x in atom_res_names
        ]
        atom_element_names = [
            _atomic_num_to_element(int(x.item())) for x in self.atom_ref_element
        ]

        pdb_atoms = []
        for atom_index in range(self.num_atoms):
            if not self.atom_exists_mask[atom_index].item():
                # skip missing atoms
                continue

            chain_tag_vocab = string.ascii_uppercase + string.ascii_lowercase
            if int(atom_asym_id[atom_index].item()) >= len(chain_tag_vocab):
                logger.warning(
                    f"Too many chains for PDB file: {atom_asym_id[atom_index].item()} -- wrapping around"
                )
            atom = PDBAtom(
                record_type="ATOM" if not self.is_ligand else "HETATM",
                atom_index=atom_index,
                atom_name=atom_names[atom_index],
                alt_loc="",
                res_name_3=atom_res_names_strs[atom_index],
                chain_tag=chain_tag_vocab[
                    int(atom_asym_id[atom_index].item()) % len(chain_tag_vocab)
                ],
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

    # @classmethod
    # def cat(cls, contexts: list["PDBContext"]) -> "PDBContext":
    #     """Concatenates multiple posebuster contexts into a single context"""
    #     cat_attrs: dict[str, Tensor] = dict()
    #     for attr in cls.__annotations__.keys():
    #         cat_attrs[attr] = torch.cat([getattr(c, attr) for c in contexts], dim=0)
    #     return cls(**cat_attrs)


def _atomic_num_to_element(atomic_num: int) -> str:
    return gemmi.Element(atomic_num).name


def entity_to_pdb_atoms(entity: PDBContext) -> list[list[PDBAtom]]:
    """Writes a single tokenized entity to PDB file"""
    pdb_atoms = entity.get_pdb_atoms()
    chains = defaultdict(list)
    for atom in pdb_atoms:
        chains[atom.chain_tag].append(atom)
    return list(chains.values())


def entities_to_pdb_file(entities: list[PDBContext], path: str):
    pdb_atoms: list[list[PDBAtom]] = []
    for entity in entities:
        pdb_atoms = pdb_atoms + entity_to_pdb_atoms(entity)
    write_pdb(pdb_atoms, path)


def pdb_context_from_batch(
    d: dict, coords: Tensor, plddt: Tensor | None = None
) -> PDBContext:
    return PDBContext(
        token_residue_index=d["token_residue_index"][0],
        token_asym_id=d["token_asym_id"][0],
        token_entity_type=d["token_entity_type"][0],
        token_residue_names=d["token_residue_name"][0],
        atom_token_index=d["atom_token_index"][0],
        atom_ref_element=d["atom_ref_element"][0],
        atom_ref_mask=d["atom_ref_mask"][0],
        atom_coords=coords[0],
        atom_exists_mask=d["atom_exists_mask"][0],
        atom_ref_name_chars=d["atom_ref_name_chars"][0],
        atom_bfactor_or_plddt=plddt[0] if plddt is not None else None,
    )


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
        for ords_atom in tensor.squeeze(0)
    ]

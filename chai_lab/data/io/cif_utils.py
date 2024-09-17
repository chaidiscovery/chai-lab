import logging
from collections import defaultdict
from pathlib import Path

import gemmi
import modelcif
import torch
from ihm import ChemComp, DNAChemComp, LPeptideChemComp, RNAChemComp
from modelcif import Assembly, AsymUnit, Entity, dumper, model
from torch import Tensor

from chai_lab.data.io import pdb_utils
from chai_lab.data.io.pdb_utils import (
    PDBContext,
    entity_to_pdb_atoms,
    get_pdb_chain_name,
    pdb_context_from_batch,
)
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import restype_3to1
from chai_lab.utils.typing import Float

logger = logging.getLogger(__name__)


class _LocalPLDDT(modelcif.qa_metric.Local, modelcif.qa_metric.PLDDT):
    name = "pLDDT"
    software = None
    description = "Predicted lddt"


def token_centre_plddts(
    context: PDBContext,
    asym_id: int,
) -> tuple[list[float], list[int]]:
    assert context.atom_bfactor_or_plddt is not None
    mask = context.token_asym_id == asym_id

    atom_idces = context.token_centre_atom_index[mask]

    # residue indices for tokens
    residue_indices = context.token_residue_index[mask]

    return context.atom_bfactor_or_plddt[atom_idces].tolist(), residue_indices.tolist()


def get_chains_metadata(context: PDBContext, entity_names) -> dict[int, AsymUnit]:
    # for each chain, get chain id, entity id, full sequence
    token_res_names = context.token_res_names_to_string

    asym_id2asym_unit = {}

    for asym_id in torch.unique(context.token_asym_id):
        if asym_id == 0:  # padding
            continue

        [token_indices] = torch.where(context.token_asym_id == asym_id)
        chain_token_res_names = [token_res_names[i] for i in token_indices]

        residue_indices = context.token_residue_index[token_indices]

        # check no missing residues
        max_res = int(torch.max(residue_indices).item())

        any_token_in_resi = residue_indices.new_zeros([max_res + 1]) - 99

        any_token_in_resi[residue_indices] = torch.arange(
            len(residue_indices),
            dtype=residue_indices.dtype,
            device=residue_indices.device,
        )

        assert any_token_in_resi.min() >= 0

        sequence = [chain_token_res_names[i] for i in any_token_in_resi]
        entity_id = context.token_entity_id[token_indices[0]]

        entity_type: int = context.get_chain_entity_type(asym_id)

        chain_id_str = get_pdb_chain_name(asym_id.item())

        asym_id2asym_unit[int(asym_id.item())] = AsymUnit(
            entity=Entity(
                # sequence is a list of ChemComponents for aminoacids/bases
                sequence=[_to_chem_component(resi, entity_type) for resi in sequence],
                description=entity_names[int(entity_id)],
            ),
            details=f"Chain {chain_id_str}",
            id=chain_id_str,
        )

    return asym_id2asym_unit


def _to_chem_component(res_name_3: str, entity_type: int):
    match entity_type:
        case EntityType.LIGAND.value:
            code = res_name_3
            return ChemComp(res_name_3, code, code_canonical=code)
        case EntityType.PROTEIN.value:
            code = restype_3to1.get(res_name_3, res_name_3)
            one_letter_code = gemmi.find_tabulated_residue(res_name_3).one_letter_code
            return LPeptideChemComp(res_name_3, code, code_canonical=one_letter_code)
        case EntityType.DNA.value:
            code = res_name_3
            # canonical code is DA -> A for DNA in cif files from wwpdb
            return DNAChemComp(res_name_3, code, code_canonical=code[-1])
        case EntityType.RNA.value:
            code = res_name_3
            return RNAChemComp(res_name_3, code, code_canonical=code)
        case _:
            raise NotImplementedError()


def context_to_cif(context: PDBContext, outpath: Path, entity_names: dict[int, str]):
    asym_id2asym_unit = get_chains_metadata(context, entity_names=entity_names)

    pdb_atoms: list[list] = entity_to_pdb_atoms(context)

    _assembly = Assembly(elements=asym_id2asym_unit.values(), name="Assembly 1")

    _model = model.AbInitioModel(_assembly, name="pred_model_1")

    for atoms_list in pdb_atoms:
        for a in atoms_list:
            _model.add_atom(
                model.Atom(
                    asym_unit=asym_id2asym_unit[a.asym_id],
                    type_symbol=a.element,
                    seq_id=int(a.residue_index),
                    atom_id=a.atom_name,
                    x=a.pos[0],
                    y=a.pos[1],
                    z=a.pos[2],
                    het=False,  # TODO should be true for ligands.
                    biso=a.b_factor,
                    occupancy=1.00,
                )
            )

    for asym_id, cif_asym_unit in asym_id2asym_unit.items():
        entity_type = context.get_chain_entity_type(asym_id)

        if entity_type != EntityType.LIGAND.value:
            # token centres are Ca or C1' for nucleic acids
            chain_plddts, residue_indices = token_centre_plddts(context, asym_id)

            for residue_idx, plddt in zip(residue_indices, chain_plddts, strict=True):
                _model.qa_metrics.append(
                    _LocalPLDDT(cif_asym_unit.residue(residue_idx + 1), plddt)
                )

    system = modelcif.System(title="Chai-1 predicted structure")
    system.authors = ["Chai Discovery team"]
    model_group = model.ModelGroup([_model], name="pred")
    system.model_groups.append(model_group)

    dumper.write(open(outpath, "w"), systems=[system])


def outputs_to_cif(
    coords: Float[Tensor, "1 n_atoms 3"],
    output_batch: dict,
    write_path: Path,
    entity_names: dict[int, str],
    bfactors: Float[Tensor, "1 n_atoms"] | None = None,
):
    context = pdb_context_from_batch(output_batch, coords, plddt=bfactors)
    write_path.parent.mkdir(parents=True, exist_ok=True)
    context_to_cif(context, write_path, entity_names)
    logger.info(f"saved cif file to {write_path}")


def outputs_to_cif_new(
    coords: Float[Tensor, "1 n_atoms 3"],
    output_batch: dict,
    write_path: Path,
    entity_names: dict[int, str],
    bfactors: Float[Tensor, "1 n_atoms"] | None = None,
):
    context = pdb_context_from_batch(output_batch, coords, plddt=bfactors)
    write_path.parent.mkdir(parents=True, exist_ok=True)
    new_context_to_cif_atoms(coords, context, entity_names, out_path=write_path)
    logger.info(f"saved cif file to {write_path}")


def new_context_to_cif_atoms(
    coords: Float[Tensor, "1 n_atoms 3"],
    context: PDBContext,
    entity_names: dict[int, str],
    out_path: Path,
):
    # TODO confirm all tensors are on CPU

    asym_id2asym_unit = get_chains_metadata(context, entity_names=entity_names)

    atom_asym_id = context.token_asym_id[context.atom_token_index]
    # atom level attributes
    atom_names = pdb_utils._tensor_to_atom_names(context.atom_ref_name_chars)

    # allows to disambiguate between atoms in logands
    asym_id_atom_name2count: dict = defaultdict(int)

    mc_assembly = Assembly(elements=asym_id2asym_unit.values(), name="Assembly 1")
    mc_model = model.AbInitioModel(mc_assembly, name="pred_model_1")

    for atom_index, token_index in enumerate(context.atom_token_index):
        if not context.atom_exists_mask[atom_index].item():
            # skip missing atoms
            continue

        x, y, z = coords[0, atom_index].tolist()
        asym_id = int(atom_asym_id[atom_index].item())
        atom_name = atom_names[atom_index]  # CA, CB, etc

        is_ligand = context.get_chain_entity_type(asym_id) == EntityType.LIGAND.value

        atom_id = atom_names[atom_index]
        if is_ligand:
            # need unique names for atoms, e.g. two carbon atoms
            asym_id_atom_name2count[asym_id, atom_name] += 1
            counter = asym_id_atom_name2count[asym_id, atom_name]
            atom_id = atom_id + f"_{counter}"

        atomic_num = context.atom_ref_element[atom_index].item()
        assert isinstance(atomic_num, int)

        mc_model.add_atom(
            model.Atom(
                asym_unit=asym_id2asym_unit[asym_id],
                type_symbol=gemmi.Element(atomic_num).name,  # e.g. H, C, O, N
                # enumerate residues from 1 in output.
                seq_id=int(context.token_residue_index[token_index].item()) + 1,
                atom_id=atom_id,
                x=x,
                y=y,
                z=z,
                het=False,  # TODO change to is_ligand
                biso=(
                    1.00
                    if context.atom_bfactor_or_plddt is None
                    else context.atom_bfactor_or_plddt[atom_index].item()
                ),
                occupancy=1.00,
            )
        )

    for asym_id, cif_asym_unit in asym_id2asym_unit.items():
        entity_type = context.get_chain_entity_type(asym_id)

        if entity_type != EntityType.LIGAND.value:
            # token centres are Ca or C1' for nucleic acids
            chain_plddts, residue_indices = token_centre_plddts(context, asym_id)

            for residue_idx, plddt in zip(residue_indices, chain_plddts, strict=True):
                mc_model.qa_metrics.append(
                    _LocalPLDDT(cif_asym_unit.residue(residue_idx + 1), plddt)
                )

    model_group = model.ModelGroup([mc_model], name="pred")

    system = modelcif.System(title="Chai-1 predicted structure")
    system.authors = ["Chai Discovery team"]
    system.model_groups.append(model_group)

    dumper.write(open(out_path, "w"), systems=[system])

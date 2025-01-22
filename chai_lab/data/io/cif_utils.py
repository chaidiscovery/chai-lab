# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
import string
from collections import defaultdict
from pathlib import Path

import gemmi
import modelcif
import torch
from einops import rearrange
from ihm import (
    DNAChemComp,
    LPeptideChemComp,
    NonPolymerChemComp,
    RNAChemComp,
    SaccharideChemComp,
)
from modelcif import Assembly, AsymUnit, Entity, dumper, model
from torch import Tensor

from chai_lab.data.io.pdb_utils import PDBContext, pdb_context_from_batch
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import restype_3to1
from chai_lab.utils.typing import Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


class _LocalPLDDT(modelcif.qa_metric.Local, modelcif.qa_metric.PLDDT):
    name = "pLDDT"
    software = None
    description = "Predicted lddt"


_CHAIN_VOCAB = [*string.ascii_uppercase, *string.ascii_lowercase]
# single-letter, then double-letter
_CHAIN_VOCAB = _CHAIN_VOCAB + [x + y for x in _CHAIN_VOCAB for y in _CHAIN_VOCAB]


def get_chain_letter(asym_id: int) -> str:
    """Get chain given a one-indexed asym_id."""
    assert asym_id > 0 and asym_id <= len(_CHAIN_VOCAB)
    vocab_index = asym_id - 1  # 1 -> A, 2 -> B
    return _CHAIN_VOCAB[vocab_index]


@typecheck
def _tensor_to_atom_names(
    tensor: Int[Tensor, "n 4"] | UInt8[Tensor, "n 4"],
) -> list[str]:
    return [
        "".join([chr(ord_val + 32) for ord_val in ords_atom]).rstrip()
        for ords_atom in tensor
    ]


@typecheck
def token_centre_plddts(
    context: PDBContext,
    plddts: Float[Tensor, "n"],
    asym_id: int,
) -> tuple[list[float], list[int]]:
    assert plddts is not None
    mask = context.token_asym_id == asym_id

    atom_idces = context.token_centre_atom_index[mask]

    # residue indices for tokens
    residue_indices = context.token_residue_index[mask]

    return plddts[atom_idces].tolist(), residue_indices.tolist()


def get_chains_metadata(
    context: PDBContext, asymid2entity_name: dict[int, str]
) -> dict[int, AsymUnit]:
    """Return mapping from asym id to AsymUnit objects."""
    assert context.asym_id2entity_type.keys() == asymid2entity_name.keys()
    # for each chain, get chain id, entity id, full sequence
    token_res_names = context.token_res_names_to_string

    asym_id2asym_unit = {}
    entity_key2ihm_entity = {}

    for asym_id, entity_type in context.asym_id2entity_type.items():
        assert asym_id != 0, "zero is padding for asym_id"

        [token_indices] = torch.where(context.token_asym_id == asym_id)
        chain_token_res_names = [token_res_names[i] for i in token_indices]

        residue_indices = context.token_residue_index[token_indices]

        # check no missing residues
        max_res = int(torch.max(residue_indices).item())

        any_token_in_resi = residue_indices.new_zeros([max_res + 1]) - 99

        any_token_in_resi[residue_indices] = torch.arange(
            len(residue_indices), dtype=any_token_in_resi.dtype
        )

        assert any_token_in_resi.min() >= 0

        sequence = [chain_token_res_names[i] for i in any_token_in_resi]

        if entity_type == EntityType.LIGAND.value:
            entity_key = (entity_type, asym_id)  # each ligand = separate entity.
        else:
            entity_key = (entity_type, *sequence)

        asym_entity_name = asymid2entity_name[asym_id]
        if entity_key not in entity_key2ihm_entity:
            entity_key2ihm_entity[entity_key] = Entity(
                # sequence is a list of ChemComponents for aminoacids/bases
                sequence=[
                    _to_chem_component(resi, entity_type, asym_id) for resi in sequence
                ],
                # will be named same as first among replicas
                description=f"Entity {asym_entity_name}",
            )

        asym_id2asym_unit[asym_id] = AsymUnit(
            entity=entity_key2ihm_entity[entity_key],
            details=f"Chain {asym_entity_name}",
            id=asym_entity_name,
        )

    return asym_id2asym_unit


def _to_chem_component(res_name_3: str, entity_type: int, asym_id: int):
    match entity_type:
        case EntityType.LIGAND.value:
            return NonPolymerChemComp(id=res_name_3 + str(asym_id))
        case EntityType.MANUAL_GLYCAN.value:
            return SaccharideChemComp(id=res_name_3, name=res_name_3)
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
            raise NotImplementedError(f"Cannot handle entity type: {entity_type}")


@typecheck
def save_to_cif(
    coords: Float[Tensor, "1 n_atoms 3"],
    output_batch: dict,
    write_path: Path,
    asym_entity_names: dict[int, str],
    bfactors: Float[Tensor, "1 n_atoms"] | None = None,
):
    write_path.parent.mkdir(parents=True, exist_ok=True)
    new_context_to_cif_atoms(
        # all inputs to function are on CPU
        coords=rearrange(coords, "1 n c -> n c", c=3).cpu(),
        plddts=None if bfactors is None else rearrange(bfactors, "1 n -> n").cpu(),
        context=pdb_context_from_batch(output_batch),
        asym_entity_names=asym_entity_names,
        out_path=write_path,
    )
    logger.info(f"saved cif file to {write_path}")


@typecheck
def new_context_to_cif_atoms(
    coords: Float[Tensor, "n_atoms 3"],
    plddts: Float[Tensor, "n_atoms"] | None,
    context: PDBContext,
    asym_entity_names: dict[int, str],
    out_path: Path,
):
    asym_id2asym_unit = get_chains_metadata(
        context, asymid2entity_name=asym_entity_names
    )

    atom_asym_id = context.token_asym_id[context.atom_token_index]
    # atom level attributes
    atom_names = _tensor_to_atom_names(context.atom_ref_name_chars)

    # allows to disambiguate between atoms in logands
    asym_id_atom_name2count: dict = defaultdict(int)

    mc_assembly = Assembly(elements=asym_id2asym_unit.values(), name="Assembly 1")
    mc_model = model.AbInitioModel(mc_assembly, name="pred_model_1")

    for atom_index, token_index in enumerate(context.atom_token_index):
        if not context.atom_exists_mask[atom_index].item():
            # skip missing atoms
            continue

        x, y, z = coords[atom_index].tolist()
        asym_id = int(atom_asym_id[atom_index].item())
        atom_name = atom_names[atom_index]  # CA, CB, etc

        is_ligand = context.asym_id2entity_type[asym_id] == EntityType.LIGAND.value

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
                het=is_ligand,
                biso=(None if plddts is None else plddts[atom_index].item()),
                occupancy=1.00,
            )
        )

    if plddts is not None:
        for asym_id, cif_asym_unit in asym_id2asym_unit.items():
            entity_type = context.asym_id2entity_type[asym_id]

            if entity_type != EntityType.LIGAND.value:
                # token centres are Ca for proteins, C1' for nucleic acids
                chain_plddts, residue_indices = token_centre_plddts(
                    context, plddts, asym_id
                )

                for residue_idx, plddt in zip(
                    residue_indices, chain_plddts, strict=True
                ):
                    mc_model.qa_metrics.append(
                        _LocalPLDDT(cif_asym_unit.residue(residue_idx + 1), plddt)
                    )

    model_group = model.ModelGroup([mc_model], name="pred")

    system = modelcif.System(title="Chai-1 predicted structure")
    system.authors = ["Chai Discovery team"]
    system.model_groups.append(model_group)

    dumper.write(open(out_path, "w"), systems=[system])

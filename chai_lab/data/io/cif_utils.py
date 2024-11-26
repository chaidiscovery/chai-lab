# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from pathlib import Path

import gemmi
import modelcif
import torch
from ihm import ChemComp, DNAChemComp, LPeptideChemComp, RNAChemComp
from modelcif import Assembly, AsymUnit, Entity, dumper, model
from torch import Tensor

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


def get_chains_metadata(context: PDBContext) -> list[dict]:
    # for each chain, get chain id, entity id, full sequence
    token_res_names = context.token_res_names_to_string

    records = []

    for asym_id in torch.unique(context.token_asym_id):
        if asym_id == 0:  # padding
            continue

        token_indices = torch.where(context.token_asym_id == asym_id)[0]
        chain_token_res_names = [token_res_names[i] for i in token_indices]

        residue_indices = context.token_residue_index[token_indices]

        # check no missing residues
        max_res = torch.max(residue_indices).item()
        tokens_in_residue = (
            residue_indices == torch.arange(max_res + 1)[..., None]
        )  # (residues, tokens)
        assert tokens_in_residue.any(dim=-1).all()

        first_token_in_resi = torch.argmax(tokens_in_residue.int(), dim=-1)

        sequence = [chain_token_res_names[i] for i in first_token_in_resi]
        entity_id = context.token_entity_id[token_indices][0]

        entity_type: int = context.get_chain_entity_type(asym_id)

        records.append(
            {
                "sequence": sequence,
                "entity_id": entity_id.item(),
                "asym_id": asym_id.item(),
                "entity_type": entity_type,
            }
        )

    return records


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
            raise NotImplementedError


def sequence_to_chem_comps(sequence: list[str], entity_type: int) -> list[ChemComp]:
    return [_to_chem_component(resi, entity_type) for resi in sequence]


def context_to_cif(context: PDBContext, outpath: Path, entity_names: dict[int, str]):
    records = get_chains_metadata(context)

    entities_map = {}
    for record in records:
        entity_id = record["entity_id"]
        if entity_id in entities_map:
            continue

        chem_components = sequence_to_chem_comps(
            record["sequence"], record["entity_type"]
        )
        cif_entity = Entity(chem_components, description=entity_names[entity_id])

        entities_map[entity_id] = cif_entity

    chains_map = {}

    def _make_chain(record: dict) -> AsymUnit:
        asym_id = record["asym_id"]
        chain_id_str = get_pdb_chain_name(asym_id)
        entity_id = record["entity_id"]

        return AsymUnit(
            entities_map[entity_id],
            details=f"Chain {chain_id_str}",
            id=chain_id_str,
        )

    chains_map = {r["asym_id"]: _make_chain(r) for r in records}

    pdb_atoms: list[list] = entity_to_pdb_atoms(context)

    _assembly = Assembly(chains_map.values(), name="Assembly 1")

    class PredModel(model.AbInitioModel):
        def get_atoms(self):
            for atoms_list in pdb_atoms:
                for a in atoms_list:
                    yield model.Atom(
                        asym_unit=chains_map[a.asym_id],
                        type_symbol=a.element,
                        seq_id=int(a.residue_index),
                        atom_id=a.atom_name,
                        x=a.pos[0],
                        y=a.pos[1],
                        z=a.pos[2],
                        het=False,
                        biso=a.b_factor,
                        occupancy=1.00,
                    )

    _model = PredModel(_assembly, name="pred_model_1")

    for asym_id, cif_asym_unit in chains_map.items():
        entity_type = context.get_chain_entity_type(asym_id)

        if entity_type != EntityType.LIGAND.value:
            # token centres are Ca or C1' for nucleic acids
            chain_plddts, residue_indices = token_centre_plddts(
                context,
                asym_id,
            )

            for residue_idx, plddt in zip(
                residue_indices,
                chain_plddts,
            ):
                _model.qa_metrics.append(
                    _LocalPLDDT(cif_asym_unit.residue(residue_idx + 1), plddt)
                )

    system = modelcif.System(title="Chai-1 predicted structure")
    system.authors = ["Chai team"]
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

import logging
import string
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import gemmi

from chai_lab.data.dataset.structure.all_atom_residue_tokenizer import (
    AllAtomResidueTokenizer,
    _make_sym_ids,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.fasta import parse_modified_fasta_sequence, read_fasta
from chai_lab.data.parsing.structure.all_atom_entity_data import AllAtomEntityData
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.parsing.structure.residue import Residue, get_restype
from chai_lab.data.residue_constants import (
    new_ligand_residue_name,
    residue_types_with_nucleotides_order,
)
from chai_lab.data.sources.rdkit import RefConformerGenerator

logger = logging.getLogger(__name__)


@dataclass
class Input:
    sequence: str
    entity_type: int


def get_lig_residues(
    smiles: str,
) -> list[Residue]:
    return [
        Residue(
            name=new_ligand_residue_name,
            label_seq=0,
            restype=residue_types_with_nucleotides_order["X"],
            residue_index=0,
            is_missing=False,
            b_factor_or_plddt=0.0,
            conformer_data=None,
            smiles=smiles,
        )
    ]


def get_polymer_residues(
    residue_names: list[str],
    entity_type: EntityType,
):
    residues = []
    for i, residue_name in enumerate(residue_names):
        residues.append(
            Residue(
                name=residue_name,
                label_seq=i,
                restype=get_restype(
                    gemmi.find_tabulated_residue(residue_name), entity_type
                ),
                residue_index=i,
                is_missing=False,
                b_factor_or_plddt=1.0,
                conformer_data=None,
            )
        )
    return residues


def _synth_subchain_id(idx: int) -> str:
    n = len(string.ascii_uppercase)
    retval = ""
    while idx >= 0:
        retval = string.ascii_uppercase[idx % n] + retval
        idx = idx // n - 1
    return retval


def raw_inputs_to_entitites_data(
    inputs: list[Input], identifier: str = "test"
) -> list[AllAtomEntityData]:
    entities = []

    # track unique entities
    entity_to_index: dict[tuple[EntityType, tuple[str, ...]], int] = {}

    for i, input in enumerate(inputs):
        # Parse residues based on entity type
        residues = None
        match entity_type := EntityType(input.entity_type):
            case EntityType.LIGAND:
                residues = get_lig_residues(smiles=input.sequence)

            case EntityType.PROTEIN | EntityType.RNA | EntityType.DNA:
                parsed_sequence: list = parse_modified_fasta_sequence(
                    input.sequence, entity_type
                )
                residues = get_polymer_residues(parsed_sequence, entity_type)
            case _:
                raise NotImplementedError
        assert residues is not None

        # Determine the entity id (unique integer for each distinct sequence)
        # NOTE very important for recognizing things like homo polymers
        seq: tuple[str, ...] = tuple(res.name for res in residues)
        entity_key: tuple[EntityType, tuple[str, ...]] = (entity_type, seq)
        if entity_key in entity_to_index:
            entity_id = entity_to_index[entity_key]
        else:
            entity_id = len(entity_to_index)
            entity_to_index[entity_key] = entity_id

        entities.append(
            AllAtomEntityData(
                residues,
                full_sequence=[residue.name for residue in residues],
                resolution=0.0,
                release_datetime=datetime.now(),
                pdb_id=identifier,
                source_pdb_chain_id=_synth_subchain_id(i),
                entity_name=f"entity_{i}_{entity_type.name}",
                entity_id=entity_id,
                method="none",
                entity_type=entity_type,
                subchain_id=_synth_subchain_id(i),
            )
        )

    assert len(entities) == len(inputs)
    return entities


def load_chains_from_raw(
    inputs: list[Input],
    identifier: str = "test",
    tokenizer: AllAtomResidueTokenizer | None = None,
) -> list[Chain]:
    """
    loads and tokenizes each input chain
    """

    if tokenizer is None:
        conformer_generator = RefConformerGenerator()
        tokenizer = AllAtomResidueTokenizer(conformer_generator)

    # Extract the entity data from the gemmi structure.
    entities: list[AllAtomEntityData] = raw_inputs_to_entitites_data(
        inputs,
        identifier=identifier,
    )

    # Tokenize the entity data
    structure_contexts = []
    sym_ids = _make_sym_ids([x.entity_id for x in entities])
    for idx, (entity_data, sym_id) in enumerate(zip(entities, sym_ids)):
        try:
            tok = tokenizer._tokenize_entity(
                entity_data,
                chain_id=idx + 1,
                sym_id=sym_id,
            )
            structure_contexts.append(tok)
        except Exception:
            logger.exception(f"Failed to tokenize input {inputs[idx]}")

    # Join the untokenized entity data with the tokenized chain data, removing
    # chains we failed to tokenize
    chains = [
        Chain(entity_data=entity_data, structure_context=structure_context)
        for entity_data, structure_context in zip(entities, structure_contexts)
        if structure_context is not None
    ]

    return chains


def read_inputs(fasta_file: str | Path, length_limit: int | None = None) -> list[Input]:
    """Read inputs from a fasta file.

    If the total length of sequences' character count is greater than length limit,
    return an empty list. Note that character count is not the same as token count, but
    is an easy approximation (smiles length is somewhat proportion to number of atoms in
    a ligand, number of residues approximates number of tokens with modified amino acids
    adding to it, etc.).
    """
    sequences = read_fasta(fasta_file)

    retval: list[Input] = []
    total_length: int = 0
    for desc, sequence in sequences:
        logger.info(f"[fasta] [{fasta_file}] {desc} {len(sequence)}")
        # get the type of the sequence
        entity_str = desc.split("|")[0].strip()
        match entity_str:
            case "protein":
                entity_type = EntityType.PROTEIN
            case "ligand":
                entity_type = EntityType.LIGAND
            case "rna":
                entity_type = EntityType.RNA
            case "dna":
                entity_type = EntityType.DNA
            case _:
                raise ValueError(f"{entity_str} is not a valid entity type")
        retval.append(Input(sequence, entity_type.value))
        total_length += len(sequence)

    if length_limit is not None and total_length > length_limit:
        logger.warning(
            f"[fasta] [{fasta_file}] too many chars ({total_length} > {length_limit}); skipping"
        )
        return []

    return retval

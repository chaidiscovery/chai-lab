# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

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
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.fasta import get_residue_name, read_fasta
from chai_lab.data.parsing.glycans import glycan_string_residues
from chai_lab.data.parsing.input_validation import (
    constituents_of_modified_fasta,
    identify_potential_entity_types,
)
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
    entity_name: str


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
) -> list[Residue]:
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
    """Load an entity for each raw input."""
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
                parsed_sequence: list | None = constituents_of_modified_fasta(
                    input.sequence
                )
                assert (
                    parsed_sequence is not None
                ), f"incorrect FASTA: {parsed_sequence=} "
                expanded_sequence = [
                    get_residue_name(r, entity_type=entity_type) if len(r) == 1 else r
                    for r in parsed_sequence
                ]
                residues = get_polymer_residues(expanded_sequence, entity_type)
            case EntityType.MANUAL_GLYCAN:
                residues = glycan_string_residues(glycan_string=input.sequence)
            case _:
                raise NotImplementedError
        assert residues is not None

        # Determine the entity id (unique integer for each distinct sequence)
        # NOTE because ligand residues have a single "LIG" residue name, the name field
        # cannot be used to distinguish them. Instead, we use the sequence field itself,
        # which should contain the SMILES string. This is not ideal, as it fails to
        # distinguish betwen different SMILES strings that represent the same molecule,
        # but should capture most cases.
        # We do not need to do special check on glycans because they are specified as a
        # string of monosaccharides, which behaves similarly to a string of amino acid
        # residues.
        seq: tuple[str, ...] = (
            (input.sequence,)
            if input.entity_type == EntityType.LIGAND.value
            else tuple(res.name for res in residues)
        )
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
                entity_name=input.entity_name,
                entity_id=entity_id,
                method="none",
                entity_type=entity_type,
                subchain_id=_synth_subchain_id(i),
                original_record=input.sequence,
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
    Loads and tokenizes each input chain; skips over inputs that fail to tokenize.
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
    structure_contexts: list[AllAtomStructureContext | None] = []
    sym_ids = _make_sym_ids([x.entity_id for x in entities])
    for entity_data, sym_id in zip(entities, sym_ids):
        # chain index should not count null contexts that result from failed tokenization
        chain_index = sum(ctx is not None for ctx in structure_contexts) + 1
        try:
            tok = tokenizer._tokenize_entity(
                entity_data,
                chain_id=chain_index,
                sym_id=sym_id,
            )
            if tok is None:
                logger.exception(f"Failed to tokenize input {entity_data=}  {sym_id=}")
        except Exception as e:
            logger.exception(
                f"Failed to tokenize input {entity_data=}  {sym_id=}", exc_info=e
            )
            tok = None
        structure_contexts.append(tok)

    # Join the untokenized entity data with the tokenized chain data, removing
    # chains we failed to tokenize
    chains = [
        Chain(entity_data=entity_data, structure_context=structure_context)
        for entity_data, structure_context in zip(
            entities, structure_contexts, strict=True
        )
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
        # examples of inputs
        # 'protein|example-of-protein'
        # 'protein|name=example-of-protein'
        # 'protein|name=example-of-protein|use_esm=true' # example how it can be in the future

        entity_str, *desc_parts = desc.split("|")

        match entity_str.lower().strip():
            case "protein":
                entity_type = EntityType.PROTEIN
            case "ligand":
                entity_type = EntityType.LIGAND
            case "rna":
                entity_type = EntityType.RNA
            case "dna":
                entity_type = EntityType.DNA
            case "glycan":
                entity_type = EntityType.MANUAL_GLYCAN
            case _:
                raise ValueError(f"{entity_str} is not a valid entity type")

        match desc_parts:
            case []:
                raise ValueError(f"label is not provided in {desc=}")
            case [label_part]:
                label_part = label_part.strip()
                if "=" in label_part:
                    field_name, entity_name = label_part.split("=")
                    assert field_name == "name"
                else:
                    entity_name = label_part
            case _:
                raise ValueError(f"Unsupported inputs: {desc=}")

        possible_types = identify_potential_entity_types(sequence)
        if len(possible_types) == 0:
            logger.error(f"Provided {sequence=} is invalid")
        elif entity_type not in possible_types:
            types_fmt = "/".join(str(et.name) for et in possible_types)
            logger.warning(
                f"Provided {sequence=} is likely {types_fmt}, not {entity_type.name}"
            )

        retval.append(Input(sequence, entity_type.value, entity_name))
        total_length += len(sequence)

    if length_limit is not None and total_length > length_limit:
        raise ValueError(
            f"[fasta] [{fasta_file}] too many chars ({total_length} > {length_limit}); skipping"
        )

    return retval

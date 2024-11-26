# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging

import gemmi

from chai_lab.data import residue_constants
from chai_lab.data.parsing.structure.entity_type import EntityType

logger = logging.getLogger(__name__)


def fasta_one_letter_sequence(residue_codes: list[str]) -> str:
    """
    Converts a list of residue names into a one-letter-code sequence
    """
    return "".join(
        [gemmi.find_tabulated_residue(res).fasta_code() for res in residue_codes]
    )


def protein_one_letter_sequence(residue_codes: list[str]) -> str:
    """
    Converts a list of protein residue names into a one-letter-code sequence.
    Probably equivalent to gemmi fasta_code() method but kept for consistency
    with old parsing + to be explicit about how non-standard res are handled (with X)
    """
    return "".join([_get_protein_only_residue_token(res) for res in residue_codes])


def protein_one_letter_sequence_with_mods(residue_codes: list[str]) -> str:
    """
    Convert a list of protein residue names into a one-letter code sequence,
    insert non-standard residues as [FOO] where FOO corresponds to the residue code of
    that non-standard residue.

    For example, 1PFH is ...APNGL[HIP]TRP... where HIP is the modified residue.
    """
    return "".join(
        [
            _get_protein_only_residue_token(res, mods_in_brackets=True)
            for res in residue_codes
        ]
    )


def _get_protein_only_residue_token(
    three_letter_code: str,
    mods_in_brackets: bool = False,
) -> str:
    """Encodes everything that is not a standard amino acid as X if nonstandard_as_X is
    True, otherwise return nonstandard FOO as [FOO]"""
    residue_info = gemmi.find_tabulated_residue(three_letter_code)
    # Standard amino acids are always given as single letters
    if residue_info.is_amino_acid() and residue_info.is_standard():
        single_letter = residue_info.one_letter_code
        single_letter = single_letter.upper()
        # non-standard residues derived from a parent std residue are lowercase
        single_letter = (
            single_letter if single_letter in residue_constants.restypes else "X"
        )
        return single_letter
    else:
        if mods_in_brackets:
            return f"[{three_letter_code}]"
        else:
            # non-standard residues derived from a parent std residue may have a
            # lowercase one-letter code; make this upper case.
            single_letter = residue_info.one_letter_code.upper()
            return single_letter if single_letter in residue_constants.restypes else "X"


def _get_residue_token(
    three_letter_code: str,
    entity_type: EntityType,
) -> str:
    """
    Encodes amino-acids and nucleic acids into corresponding tokens
    20 standard AAs + X
    4 RNA bases + RX
    4 DNA bases + DX
    """
    residue_info = gemmi.find_tabulated_residue(three_letter_code)
    if residue_info.is_amino_acid():
        single_letter = residue_info.one_letter_code
        single_letter = single_letter.upper()
        # non-standard residues derived from a parent std residue are lowercase
        single_letter = (
            single_letter if single_letter in residue_constants.restypes else "X"
        )
        return single_letter

    elif residue_info.is_nucleic_acid() and entity_type == EntityType.RNA:
        return "R{}".format(residue_info.one_letter_code)

    elif residue_info.is_nucleic_acid() and entity_type == EntityType.DNA:
        return "D{}".format(residue_info.one_letter_code)

    else:
        # more properties at https://gemmi.readthedocs.io/en/latest/mol.html#built-in-data
        return "X"


def get_residue_codes(subchain: gemmi.ResidueSpan, entity: gemmi.Entity) -> list[str]:
    """
    Get list of residue codes (3-letter for protein residues,
    1 to 3 letters/digits for ligands, 1 or 2 letters for RNA/DNA)
    for a gemmi subchain
    """
    # entity.full_sequence comes from SEQRES, so it might be missing in PDB files
    if entity.full_sequence is not None and len(entity.full_sequence) > 0:
        return [
            gemmi.Entity.first_mon(item)  # Ignore point mutations
            for item in entity.full_sequence
        ]
    # this infers the sequence from the set of residues in the structure
    return [res.name for res in subchain.first_conformer()]


def tokenize_sequence(
    subchain: gemmi.ResidueSpan, entity: gemmi.Entity, entity_type: EntityType
) -> list[str]:
    three_letter_sequence = get_residue_codes(subchain, entity)

    match entity_type:
        case EntityType.PROTEIN:
            return [
                _get_protein_only_residue_token(three_letter_code)
                for three_letter_code in three_letter_sequence
            ]
        case EntityType.RNA | EntityType.DNA:
            return [
                _get_residue_token(three_letter_code, entity_type)
                for three_letter_code in three_letter_sequence
            ]
        case _:
            raise NotImplementedError

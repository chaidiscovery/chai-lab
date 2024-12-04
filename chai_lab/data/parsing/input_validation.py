# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

"""
Simple heuristics that can help with identification of EntityType
"""

import string
from string import ascii_letters

from chai_lab.data.parsing.structure.entity_type import EntityType


def constituents_of_modified_fasta(x: str) -> list[str] | None:
    """
    Accepts RNA/DNA inputs: 'agtc', 'AGT(ASP)TG', etc. Does not accept SMILES strings.
    Returns constituents, e.g, [A, G, T, ASP, T, G] or None if string is incorrect.
    Everything in returned list is single character, except for blocks specified in brackets.
    """
    x = x.strip().upper()
    # it is a bit strange that digits are here, but [NH2] was in one protein
    allowed_chars = ascii_letters + "()" + string.digits
    if not all(letter in allowed_chars for letter in x):
        return None

    current_modified: str | None = None

    constituents = []
    for letter in x:
        if letter == "(":
            if current_modified is not None:
                return None  # double open bracket
            current_modified = ""
        elif letter == ")":
            if current_modified is None:
                return None  # closed without opening
            if len(current_modified) <= 1:
                return None  # empty modification: () or single (K)
            constituents.append(current_modified)
            current_modified = None
        else:
            if current_modified is not None:
                current_modified += letter
            else:
                if letter not in ascii_letters:
                    return None  # strange single-letter residue
                constituents.append(letter)
    if current_modified is not None:
        return None  # did not close bracket
    return constituents


def identify_potential_entity_types(sequence: str) -> list[EntityType]:
    """
    Provided FASTA sequence or smiles, lists which entities those could be.
    Returns an empty list if sequence is invalid for all entity types.
    """
    sequence = sequence.strip()
    if len(sequence) == 0:
        return []
    possible_entity_types = []

    constituents = constituents_of_modified_fasta(sequence)
    if constituents is not None:
        # this can be RNA/DNA/protein.
        one_letter_constituents = set(x for x in constituents if len(x) == 1)
        if set.issubset(one_letter_constituents, set("AGTC")):
            possible_entity_types.append(EntityType.DNA)
        if set.issubset(one_letter_constituents, set("AGUC")):
            possible_entity_types.append(EntityType.RNA)
        if "U" not in one_letter_constituents:
            possible_entity_types.append(EntityType.PROTEIN)

    ascii_symbols = string.ascii_letters + string.digits + ".-+=#$%:/\\[]()<>@"
    if set.issubset(set(sequence.upper()), set(ascii_symbols)):
        possible_entity_types.extend([EntityType.LIGAND, EntityType.MANUAL_GLYCAN])
    return possible_entity_types

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from pathlib import Path
from typing import NamedTuple

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import restype_1to3_with_x

logger = logging.getLogger(__name__)

Fasta = NamedTuple("Fasta", [("header", str), ("sequence", str)])


nucleic_acid_1_to_name: dict[tuple[str, EntityType], str] = {
    ("A", EntityType.RNA): "A",
    ("U", EntityType.RNA): "U",
    ("G", EntityType.RNA): "G",
    ("C", EntityType.RNA): "C",
    ("A", EntityType.DNA): "DA",
    ("T", EntityType.DNA): "DT",
    ("G", EntityType.DNA): "DG",
    ("C", EntityType.DNA): "DC",
}


def read_fasta(file_path: str | Path) -> list[Fasta]:
    from Bio import SeqIO

    fasta_sequences = SeqIO.parse(open(file_path), "fasta")
    return [Fasta(fasta.description, str(fasta.seq)) for fasta in fasta_sequences]


def get_residue_name(
    fasta_code: str,
    entity_type: EntityType,
) -> str:
    if len(fasta_code) != 1:
        raise ValueError("Cannot handle non-single chars: {}".format(fasta_code))
    match entity_type:
        case EntityType.PROTEIN:
            return restype_1to3_with_x.get(fasta_code, "UNK")
        case EntityType.RNA | EntityType.DNA:
            # under nucleic_acid_1_to_name, DNA is mapped to D_ and RNA to _
            unk = "X" if entity_type == EntityType.RNA else "DX"
            return nucleic_acid_1_to_name.get((fasta_code, entity_type), unk)
        case _:
            raise ValueError(f"Invalid polymer entity type {entity_type}")

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from io import StringIO
from pathlib import Path
from typing import NamedTuple, Sequence, TextIO

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


def fastas_to_str(fastas: Sequence[Fasta]) -> str:
    return "".join(f">{fasta.header}\n{fasta.sequence}\n" for fasta in fastas)


def read_fasta(file_path: str | Path) -> list[Fasta]:
    with open(file_path) as source:
        return read_fasta_content(source)


def read_fasta_content(content: StringIO | TextIO) -> list[Fasta]:
    from Bio import SeqIO

    fasta_sequences = SeqIO.parse(content, "fasta")
    return [Fasta(fasta.description, str(fasta.seq)) for fasta in fasta_sequences]


def read_fasta_unique(p: Path) -> tuple[list[str], list[bytes]]:
    """Read unique sequences from the given fasta file."""
    known = set()
    descriptions: list[str] = []
    aligned_seqs: list[bytes] = []
    for rec_id, rec_seq in read_fasta(p):
        aligned_seq = rec_seq.encode()
        if aligned_seq in known:
            continue  # skip duplicate sequences
        known.add(aligned_seq)

        descriptions.append(rec_id)
        aligned_seqs.append(aligned_seq)
    return descriptions, aligned_seqs


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


def write_fastas(
    fastas: Sequence[Fasta],
    output_path: str,
):
    logger.debug(f"Writing {len(fastas)} sequences to {output_path}")
    with open(output_path, "w") as fp:
        fp.write(fastas_to_str(fastas))

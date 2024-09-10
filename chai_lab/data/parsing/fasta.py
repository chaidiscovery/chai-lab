import logging
import re
from pathlib import Path
from typing import Iterable

from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.residue_constants import restype_1to3_with_x

logger = logging.getLogger(__name__)

Fasta = tuple[str, str]
Fastas = list[Fasta]


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


def read_fasta(file_path: str | Path) -> Iterable[Fasta]:
    from Bio import SeqIO

    fasta_sequences = SeqIO.parse(open(file_path), "fasta")
    return [(fasta.id, str(fasta.seq)) for fasta in fasta_sequences]


def get_residue_name(
    fasta_code: str,
    entity_type: EntityType,
) -> str:
    match entity_type:
        case EntityType.PROTEIN:
            return restype_1to3_with_x.get(fasta_code, "UNK")
        case EntityType.RNA | EntityType.DNA:
            # under nucleic_acid_1_to_name, DNA is mapped to D_ and RNA to _
            unk = "X" if entity_type == EntityType.RNA else "DX"
            return nucleic_acid_1_to_name.get((fasta_code, entity_type), unk)
        case _:
            raise ValueError(f"Invalid polymer entity type {entity_type}")


def parse_modified_fasta_sequence(sequence: str, entity_type: EntityType) -> list[str]:
    """
    Parses a fasta-like string containing modified residues in
    brackets, returns a list of residue codes
    """
    pattern = r"[A-Z]|\[[A-Z0-9]+\]"
    residues = re.findall(pattern, sequence)

    # get full residue name if regular fasta code (not in brackets),
    # otherwise return what user passed in brackets
    parsed_residues = [
        get_residue_name(x, entity_type) if not x.startswith("[") else x.strip("[]")
        for x in residues
    ]
    return parsed_residues

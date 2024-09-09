import logging
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property

from chai_lab.data.parsing.structure import sequence
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.parsing.structure.residue import Residue
from chai_lab.data.residue_constants import standard_residue_pdb_codes
from chai_lab.utils.typing import typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class AllAtomEntityData:
    residues: list[Residue]
    full_sequence: list[str]
    resolution: float
    release_datetime: datetime | None  # None if no date found
    pdb_id: str
    source_pdb_chain_id: str
    # Unique string identifying the entity.
    entity_name: str
    # Unique integer identifying the entity, starting at 0. There is a 1:1 mapping
    # between entity_name and entity_index.
    entity_id: int
    method: str
    entity_type: EntityType
    subchain_id: str
    is_d_polypeptide: bool = False  # NOTE (mostly) exists for eval set construction

    def __post_init__(self):
        assert (
            len(self.residues) == len(self.full_sequence)
        ), f"{self.__class__.__name__} residues and full_sequence must be the same length"

    @property
    def missing_residues(self) -> list[Residue]:
        """
        Returns a list of missing residues in the entity
        """
        return [residue for residue in self.residues if residue.is_missing]

    @cached_property
    def has_modifications(self) -> bool:
        """
        Returns True if the entity has modifications; this only applies to polymers so
        is always False for ligands, waters, and unknowns.
        """
        if self.entity_type not in (
            EntityType.PROTEIN,
            EntityType.RNA,
            EntityType.DNA,
            EntityType.POLYMER_HYBRID,
        ):
            return False

        return any(res.name not in standard_residue_pdb_codes for res in self.residues)

    @property
    def is_distillation(self) -> bool:
        return self.pdb_id.startswith("AF-")

    @property
    def sequence(self) -> str:
        """Sequence with modified residues encoded as X."""
        return sequence.protein_one_letter_sequence(self.full_sequence)

    @property
    def sequence_with_mods(self) -> str:
        """Sequence with modifications encoded as [FOO] where FOO is modified residue."""
        return sequence.protein_one_letter_sequence_with_mods(self.full_sequence)

    def __str__(self) -> str:
        fields = ", ".join(
            [
                f"pdb_id={self.pdb_id}",
                f"source_pdb_chain_id={self.source_pdb_chain_id}",
                f"entity_name={self.entity_name}",
                f"entity_id={self.entity_id}",
                f"entity_type={self.entity_type}",
                f"subchain_id={self.subchain_id}",
            ]
        )
        return f"AllAtomEntityData({fields})"

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import itertools
import logging
from dataclasses import dataclass
from datetime import datetime
from functools import cached_property
from typing import Iterator

import gemmi
from typing_extensions import assert_never

from chai_lab.data.parsing.structure import sequence
from chai_lab.data.parsing.structure.entity_type import EntityType, get_entity_type
from chai_lab.data.parsing.structure.residue import Residue, get_residues
from chai_lab.data.residue_constants import standard_residue_pdb_codes
from chai_lab.utils.typing import typecheck

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class StructureMetadata:
    resolution: float
    release_date: str
    pdb_id: str
    method: str


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
    original_record: str = ""  # NOTE for glycan parsing

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


def _safe_get(info: gemmi.InfoMap, key: str, default: str = "") -> str:
    if key in info:
        return info[key]
    return default


def _is_unknown_ligand(subchain: gemmi.ResidueSpan) -> bool:
    return len(subchain) == 1 and subchain[0].name == "UNL"


def get_structure_metadata(structure: gemmi.Structure) -> StructureMetadata:
    return StructureMetadata(
        resolution=structure.resolution,
        pdb_id=_safe_get(structure.info, "_entry.id"),
        release_date=_safe_get(
            structure.info, "_pdbx_database_status.recvd_initial_deposition_date"
        ),
        method=_safe_get(structure.info, "_exptl.method"),
    )


def process_entity(
    subchain: gemmi.ResidueSpan,
    structure: gemmi.Structure,
    entity_to_index: dict[str, int],
    source_pdb_chain_id: str = "A",
) -> AllAtomEntityData | None:
    """
    Process a polymer entity: protein, RNA, DNA, or hybrid, as well as ligands
    """
    gemmi_entity = structure.get_entity_of(subchain)
    entity_type = get_entity_type(gemmi_entity)
    structure_metadata = get_structure_metadata(structure)

    # used to add missing residues
    full_sequence: list[str] = sequence.get_residue_codes(subchain, gemmi_entity)

    if entity_type == EntityType.LIGAND:
        for idx, residue in enumerate(subchain):
            residue.label_seq = idx + 1

    match entity_type:
        case (
            EntityType.PROTEIN
            | EntityType.RNA
            | EntityType.DNA
            | EntityType.POLYMER_HYBRID
            | EntityType.LIGAND
        ):
            if not all(r.label_seq is not None for r in subchain):
                pdb_id = structure_metadata.pdb_id
                # All entities should have label seq at this stage
                # since we assigned it manually for ligands earlier
                logger.error(
                    f"{pdb_id=}, {source_pdb_chain_id=}, gemmi subchain {subchain.subchain_id()}"
                    "has missing label seq, skipping"
                )
                return None

            residues = get_residues(
                subchain=subchain,
                full_sequence=full_sequence,
                entity_type=entity_type,
            )

            # unique entity_id for unknown ligands
            entity_name = (
                f"UNL_{subchain.subchain_id()}"
                if _is_unknown_ligand(subchain)
                else gemmi_entity.name
            )

            if entity_name in entity_to_index:
                entity_id = entity_to_index[entity_name]
            else:
                entity_id = len(entity_to_index)
                entity_to_index[entity_name] = entity_id

            return AllAtomEntityData(
                residues=residues,
                full_sequence=full_sequence,
                entity_name=entity_name,
                entity_id=entity_id,
                entity_type=entity_type,
                source_pdb_chain_id=source_pdb_chain_id,
                subchain_id=subchain.subchain_id(),
                resolution=structure_metadata.resolution,
                release_datetime=None,
                pdb_id=structure_metadata.pdb_id,
                method=structure_metadata.method,
                is_d_polypeptide=gemmi_entity == gemmi.PolymerType.PeptideD,
            )

        case EntityType.MANUAL_GLYCAN:
            raise ValueError("Should never have manual glycans for Gemmi residues")
        case EntityType.WATER | EntityType.UNKNOWN:
            # skip waters for now
            return None

    assert_never(entity_type)


def process_chain(
    chain: gemmi.Chain,
    structure: gemmi.Structure,
    entity_to_index: dict[str, int],
) -> Iterator[AllAtomEntityData]:
    """
    Iterate over all (potentially >1) chemical entities contained in a PDB chain
    and load them into their own EntityData object
    """

    for subchain in chain.subchains():
        parsed_entity = process_entity(
            subchain,
            structure,
            entity_to_index,
            source_pdb_chain_id=chain.name,
        )
        if parsed_entity is not None:
            yield parsed_entity


def parse_structure(
    structure: gemmi.Structure,
    chain_ids: list[str] | None = None,
    make_assembly: bool = False,
) -> Iterator[Iterator[AllAtomEntityData]]:
    """
    Returns an iterator for each chain, with an iterator for each entity in the chain.
    """
    if structure.input_format == gemmi.CoorFormat.Pdb:
        raise NotImplementedError
    structure.remove_waters()

    if make_assembly:
        # make assembly involves duplicating chains in the case of homomers
        # and renaming them
        how = gemmi.HowToNameCopiedChain.Dup
        try:
            structure.transform_to_assembly(assembly_name="1", how=how)
        except RuntimeError:
            logger.warning(
                "Make assembly failed, likely no assemblies found in structure"
            )

    if chain_ids is not None:
        chains = [chain for chain in structure[0] if chain.name in chain_ids]
    else:
        chains = [chain for chain in structure[0]]

    entities_to_index: dict[str, int] = {}

    for chain in chains:
        yield process_chain(
            chain,
            structure,
            entities_to_index,
        )


def structure_to_entities_data(
    structure: gemmi.Structure,
    chain_ids: list[str] | None = None,
    subchain_ids: list[str] | None = None,
    make_assembly: bool = False,
) -> list[AllAtomEntityData]:
    entities_data = parse_structure(
        structure=structure,
        chain_ids=chain_ids,
        make_assembly=make_assembly,
    )

    if subchain_ids:
        sids = set(subchain_ids)
        entities_data = (
            (e for e in entities if e.subchain_id in sids) for entities in entities_data
        )

    results = list(itertools.chain.from_iterable(entities_data))
    return results

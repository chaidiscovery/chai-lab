# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import gemmi

from chai_lab.data.parsing.structure.all_atom_entity_data import (
    structure_to_entities_data,
)
from chai_lab.data.parsing.structure.entity_type import EntityType


def _minimal_peptide_cif(entry_id: str, poly_type: str, res1: str, res2: str) -> str:
    return f"""data_{entry_id}
_entry.id {entry_id}
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
1 {poly_type} A
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 {res1}
1 2 {res2}
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N  N  . {res1} A 1 1 ? 0.000 0.000 0.000 1.00 20.00 1 {res1} A N  1
ATOM 2 C  CA . {res1} A 1 1 ? 1.458 0.000 0.000 1.00 20.00 1 {res1} A CA 1
ATOM 3 C  C  . {res1} A 1 1 ? 2.009 1.420 0.000 1.00 20.00 1 {res1} A C  1
ATOM 4 O  O  . {res1} A 1 1 ? 1.251 2.390 0.000 1.00 20.00 1 {res1} A O  1
ATOM 5 N  N  . {res2} A 1 2 ? 3.320 1.520 0.000 1.00 20.00 2 {res2} A N  1
ATOM 6 C  CA . {res2} A 1 2 ? 4.000 2.810 0.000 1.00 20.00 2 {res2} A CA 1
ATOM 7 C  C  . {res2} A 1 2 ? 5.500 2.650 0.000 1.00 20.00 2 {res2} A C  1
ATOM 8 O  O  . {res2} A 1 2 ? 6.020 1.550 0.000 1.00 20.00 2 {res2} A O  1
"""


def test_is_d_polypeptide_detected_from_polymer_type() -> None:
    d_doc = gemmi.cif.read_string(
        _minimal_peptide_cif("TESTD", "polypeptide(D)", "DAL", "DVA")
    )
    d_entities = structure_to_entities_data(
        gemmi.make_structure_from_block(d_doc.sole_block())
    )
    assert len(d_entities) == 1
    assert d_entities[0].entity_type == EntityType.PROTEIN
    assert d_entities[0].is_d_polypeptide is True

    l_doc = gemmi.cif.read_string(
        _minimal_peptide_cif("TESTL", "polypeptide(L)", "ALA", "VAL")
    )
    l_entities = structure_to_entities_data(
        gemmi.make_structure_from_block(l_doc.sole_block())
    )
    assert len(l_entities) == 1
    assert l_entities[0].entity_type == EntityType.PROTEIN
    assert l_entities[0].is_d_polypeptide is False

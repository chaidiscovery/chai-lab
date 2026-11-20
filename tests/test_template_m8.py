# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from pathlib import Path

import gemmi
import pytest

from chai_lab.data.dataset.inference_dataset import read_inputs
from chai_lab.data.parsing.templates.m8 import parse_m8_file, parse_m8_to_template_hits
from chai_lab.tools.kalign import KalignAlignment


@pytest.mark.parametrize("query_name", ["001", "NA", "null", "protein"])
def test_template_lookup_preserves_fasta_names(
    tmp_path: Path, monkeypatch, query_name: str
):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(f">protein|name={query_name}\nAG\n")
    [chain] = read_inputs(fasta)
    m8 = tmp_path / "templates.m8"
    m8.write_text(f"{query_name}\t1abc_A\t100\t2\t0\t0\t1\t2\t1\t2\t1e-10\t20\n")

    # Local structure with a two-residue polymer; no RCSB or model call is needed.
    structure = gemmi.Structure()
    structure.name = "test"
    model = gemmi.Model(1)
    template_chain = gemmi.Chain("A")
    for i, name in enumerate(["ALA", "GLY"], start=1):
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(i, " ")
        residue.entity_type = gemmi.EntityType.Polymer
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(3.8 * i, 0, 0)
        residue.add_atom(atom)
        template_chain.add_residue(residue)
    model.add_chain(template_chain)
    structure.add_model(model)
    structure.setup_entities()
    cif = tmp_path / "template.cif"
    structure.make_mmcif_document().write_file(str(cif))

    monkeypatch.setattr(
        "chai_lab.data.parsing.templates.m8.download_cif_file",
        lambda *args, **kwargs: cif,
    )
    monkeypatch.setattr(
        "chai_lab.data.parsing.templates.m8.kalign_query_to_reference",
        lambda **kwargs: KalignAlignment("AG", "AG"),
    )

    hits = list(parse_m8_to_template_hits(chain.entity_name, chain.sequence, m8))

    assert len(hits) == 1
    assert hits[0].query_pdb_id == query_name
    assert hits[0].hit_sequence == "AG"
    assert hits[0].query_start_end == (0, 1)
    assert parse_m8_file(m8)["query_id"].tolist() == [query_name]


def test_empty_template_table(tmp_path: Path):
    m8 = tmp_path / "empty.m8"
    m8.touch()
    assert parse_m8_file(m8).empty
    assert list(parse_m8_to_template_hits("protein", "AG", m8)) == []

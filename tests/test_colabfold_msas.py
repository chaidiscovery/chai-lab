# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from pathlib import Path
from tempfile import TemporaryDirectory

from chai_lab.data.dataset.msas.colabfold import generate_colabfold_msas
from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.parsing.msas.aligned_pqt import (
    hash_sequence,
    parse_aligned_pqt_to_msa_context,
)


def test_monomer():
    sequence = "GPKAMAETRECIYYNANWELERTNQSGLERCEGEQDKRLHCYASWRNSSGTIELVKKGCWLDDFNCYDRQECVATEENPQVYFCCCEGNFCNERFTHLP"

    with TemporaryDirectory() as tmpdir:
        seq_to_aligned_pqt = generate_colabfold_msas(
            [sequence],
            msa_dir=Path(tmpdir),
            msa_server_url="https://api.colabfold.com",
            search_templates=False,
        )
        assert len(seq_to_aligned_pqt) == 1
        assert next(iter(seq_to_aligned_pqt.keys())) == sequence

        (pqt_file,) = seq_to_aligned_pqt.values()
        assert pqt_file.exists()
        ctx: MSAContext = parse_aligned_pqt_to_msa_context(pqt_file)

        # First row of the context should correspond to
        assert ctx.ith_sequence(0) == sequence


def test_multimer():
    sequences = [
        "GPKAMAETRECIYYNANWELERTNQSGLERCEGEQDKRLHCYASWRNSSGTIELVKKGCWLDDFNCYDRQECVATEENPQVYFCCCEGNFCNERFTHLP",
        "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSSYINWVRQAPGQGLEWMGTINPVSGSTSYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGGWFDYWGQGTLVTVSSHHHHHH",
    ]

    with TemporaryDirectory() as tmpdir:
        seq_to_aligned_pqt = generate_colabfold_msas(
            sequences,
            msa_dir=Path(tmpdir),
            msa_server_url="https://api.colabfold.com",
            search_templates=False,
        )
        assert len(seq_to_aligned_pqt) == 2
        for seq in sequences:
            pqt_file = seq_to_aligned_pqt[seq]
            ctx: MSAContext = parse_aligned_pqt_to_msa_context(pqt_file)
            assert ctx.depth > 1
            assert ctx.ith_sequence(0) == seq


def test_template_query_ids_are_remapped_offline(tmp_path: Path, monkeypatch):
    sequences = ["AG", "GA"]
    template_hits = tmp_path / "pdb70.m8"
    template_hits.write_text(
        "101\t1abc_A\t100\t2\t0\t0\t1\t2\t1\t2\t1e-10\t20\n"
        "102\t2def_B\t100\t2\t0\t0\t1\t2\t1\t2\t1e-10\t20\n"
    )

    def mock_mmseqs2(protein_seqs, prefix, *, use_pairing, use_templates, **kwargs):
        assert protein_seqs == sequences
        if use_pairing:
            assert not use_templates
            return ["", ""], None
        assert use_templates
        return [f">{i}\n{seq}\n" for i, seq in enumerate(sequences, 101)], str(
            template_hits
        )

    def deny_network(*args, **kwargs):
        raise AssertionError("This test must run offline")

    monkeypatch.setattr("socket.socket.connect", deny_network)
    monkeypatch.setattr(
        "chai_lab.data.dataset.msas.colabfold._run_mmseqs2", mock_mmseqs2
    )
    msa_dir = tmp_path / "msas"
    msa_dir.mkdir()
    msa_paths = generate_colabfold_msas(
        sequences,
        msa_dir=msa_dir,
        msa_server_url="https://unused.invalid",
        search_templates=True,
    )

    rows = [
        line.split("\t")
        for line in (msa_dir / "all_chain_templates.m8").read_text().splitlines()
    ]
    assert [(row[0], row[1]) for row in rows] == [
        (hash_sequence("AG"), "1abc_A"),
        (hash_sequence("GA"), "2def_B"),
    ]
    for sequence in sequences:
        ctx = parse_aligned_pqt_to_msa_context(msa_paths[sequence])
        assert ctx.ith_sequence(0) == sequence

# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from pathlib import Path
from tempfile import TemporaryDirectory

from chai_lab.data.dataset.msas.colabfold import generate_colabfold_msas
from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.parsing.msas.aligned_pqt import parse_aligned_pqt_to_msa_context


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

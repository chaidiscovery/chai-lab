# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from pathlib import Path

from chai_lab.chai1 import make_all_atom_feature_context


def _atom_label(context, atom_index: int) -> tuple[int, str]:
    token_index = context.atom_token_index[atom_index].item()
    residue_index = context.token_residue_index[token_index].item()
    return residue_index, context.atom_ref_name[atom_index]


def test_inline_modified_peptide_gets_covalent_bonds(tmp_path: Path):
    fasta_file = tmp_path / "input.fasta"
    fasta_file.write_text(">protein|name=peptide\n(TPO)(7RX)P\n")

    feature_context = make_all_atom_feature_context(
        fasta_file,
        output_dir=tmp_path / "output",
        use_esm_embeddings=False,
    )
    context = feature_context.structure_context
    left, right = context.atom_covalent_bond_indices
    bonds = {
        (_atom_label(context, atom_a.item()), _atom_label(context, atom_b.item()))
        for atom_a, atom_b in zip(left, right, strict=True)
    }

    assert bonds == {
        ((0, "C"), (1, "N08")),
        ((1, "N08"), (0, "C")),
        ((1, "C02"), (2, "N")),
        ((2, "N"), (1, "C02")),
    }

    oxy = [
        atom_index
        for atom_index in range(context.num_atoms)
        if _atom_label(context, atom_index) == (1, "OXY")
    ]
    assert len(oxy) == 1
    assert not context.atom_exists_mask[oxy[0]]


def test_terminal_modified_residue_keeps_unused_leaving_atom(tmp_path: Path):
    fasta_file = tmp_path / "input.fasta"
    fasta_file.write_text(">protein|name=peptide\nP(7RX)\n")

    feature_context = make_all_atom_feature_context(
        fasta_file,
        output_dir=tmp_path / "output",
        use_esm_embeddings=False,
    )
    context = feature_context.structure_context
    oxy = [
        atom_index
        for atom_index in range(context.num_atoms)
        if _atom_label(context, atom_index) == (1, "OXY")
    ]

    assert len(oxy) == 1
    assert context.atom_exists_mask[oxy[0]]


def test_non_peptide_ccd_is_not_automatically_linked(tmp_path: Path):
    fasta_file = tmp_path / "input.fasta"
    fasta_file.write_text(">protein|name=protein\nA(ATP)A\n")

    feature_context = make_all_atom_feature_context(
        fasta_file,
        output_dir=tmp_path / "output",
        use_esm_embeddings=False,
    )
    left, right = feature_context.structure_context.atom_covalent_bond_indices

    assert left.numel() == right.numel() == 0

from chai_lab.data.parsing.input_validation import (
    constituents_of_modified_fasta,
    identify_potential_entity_types,
)
from chai_lab.data.parsing.structure.entity_type import EntityType

from .example_inputs import example_dna, example_ligangs, example_proteins, example_rna


def test_debug():
    constituents_of_modified_fasta("[KCJ][SEP][PPN][B3S][BAL][PPN]K[NH2]")


def test_parsing():
    for ligand in example_ligangs:
        assert EntityType.LIGAND in identify_potential_entity_types(ligand)

    for protein in example_proteins:
        assert EntityType.PROTEIN in identify_potential_entity_types(protein)

    for dna in example_dna:
        assert EntityType.DNA in identify_potential_entity_types(dna)

    for rna in example_rna:
        assert EntityType.RNA in identify_potential_entity_types(rna)

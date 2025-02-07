# Copyright (c) 2024 Chai Discovery, Inc.

# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

example_ligands = [
    "C",
    "O",
    "C(C1C(C(C(C(O1)O)O)O)O)O",
    "[O-]S(=O)(=O)[O-]",
    "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C",
    "CCC1=C(c2cc3c(c(c4n3[Mg]56[n+]2c1cc7n5c8c(c9[n+]6c(c4)C(C9CCC(=O)OC/C=C(\C)/CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C)[C@H](C(=O)c8c7C)C(=O)OC)C)C=C)C=O",
    r"C=CC1=C(C)/C2=C/c3c(C)c(CCC(=O)O)c4n3[Fe@TB16]35<-N2=C1/C=c1/c(C)c(C=C)/c(n13)=C/C1=N->5/C(=C\4)C(CCC(=O)O)=C1C",
    # different ions
    "[Mg+2]",
    "[Na+]",
    "[Cl-]",
]

example_proteins = [
    "AGSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVR",
    "(KCJ)(SEP)(PPN)(B3S)(BAL)(PPN)K(NH2)",
    "XDHPX",
]


example_rna = [
    "AGUGGCUA",
    "AAAAAA",
    "AGUC",
]

example_dna = [
    "AGTGGCTA",
    "AAAAAA",
    "AGTC",
]

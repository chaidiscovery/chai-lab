# Working with bond restraints

Chai-1 supports specifying covalent bonds as input, which specify covalent linkages between atoms in the folded complex. This is useful for specifying covalent modifications such as glycosylation events, which we demonstrate below, but can be generally used to specify arbitrary "non-canonical" bonds in a structure. 

A few notes:
- Chai-1 was not trained on disulfide bonds, and we have not evaluated whether specifying such bond information yields expected behaviors. 
- These bond restraints should not be used to specify modified amino acids that already have an associated CCD code; for these examples, include the modified residue's CCD code in parentheses directly in the sequence in place of its canonical residue, e.g., `RKDES(MSE)EES` to specify a selenomethionine at the 6th position.

## Glycans

We adopt an abbreviated syntax for specifying glycans, which is best explained with a series of examples.

### Single-ring glycan

Let's say we have a glycan that is a single sugar ring, a 2-acetamido-2-deoxy-beta-D-glucopyranose. The [CCD code](https://www.rcsb.org/ligand/NAG) for this sugar is `NAG`, so we simply specify this sugar with the following fasta entry:
```
>protein|example-protein
...N...
>glycan|example-single-sugar
NAG
```

Now, a glycan is also covalently bound to a residue; to specify this, we include the following line in our restraints file (see our documentation on restraints as well):

chainA|res_idxA|chainB|res_idxB|connection_type|confidence|min_distance_angstrom|max_distance_angstrom|comment|restraint_id
|---|---|---|---|---|---|---|---|---|---|
A|N436@N|B|@C4|covalent|1.0|0.0|0.0|protein-glycan|bond1

Breaking this down, this specifies that the within chain A (the first entry in the fasta), the "N" residue at the 436-th position (1-indexed) as indicated by the "N436" prefix is bound, via its nitrogen "N" atom as indicated by the "@N" suffix, to the C4 atom in the first glycan ("@C4"). Ring numbering follows standard glycan ring number schemas. For other ligands, you will need check how the specific version of `rdkit` that we use in `chai-lab` (run `uv pip list | grep rdkit` for version) assigns atom names and use the same atom names to specify your bonds. In addition, note that the min and max distance fields are ignored, as is the confidence field. 


### Multi-ring glycan

Working through a more complex example, let's say we have a two-ring ligand such as that shown in the PDB structure [1AC5](https://www.rcsb.org/structure/1ac5). We introduce syntax for specifying bonds **within glycans** within the fasta record as such:

```
>protein|example-protein
...N...
>glycan|example-dual-sugar
NAG(1-4 NAG)
```

This syntax specifies that the root of the glycan is the leading `NAG` ring. The parentheses indicate that we are attaching another molecule to the ring directly preceding the parentheses. The `1-4` syntax "draws" a bond between the C1 atom of the previous "root" `NAG` and the C4 atom of the subsequent `NAG` ring. To specify how this glycan ought to be connected to the protein, we again use the restraints file to specify a residue and atom to which the glycan is bound, and the carbon atom within the root glycan ring that is bound.

chainA|res_idxA|chainB|res_idxB|connection_type|confidence|min_distance_angstrom|max_distance_angstrom|comment|restraint_id
|---|---|---|---|---|---|---|---|---|---|
A|N436@N|B|@C4|covalent|1.0|0.0|0.0|protein-glycan|bond1

You can chain this syntax to create longer ligands:
```
>glycan|4-NAG-in-a-linear-chain
NAG(1-4 NAG(1-4 NAG(1-4 NAG)))
```

...and to create branched ligands
```
>glycan|branched-glycan
NAG(1-4 NAG(1-4 NAG))(3-4 MAN)
```
This branched example has a root `NAG` ring with a branch with two more `NAG` rings and a branch with a single `MAN` ring. For additional examples, please refer to the examples tested in the `tests/test_glycans.py` test file.

### Example

We have included an example of how glycans can be specified under `predict_glycosylated.py` in this directory, along with its corresponding `bonds.restraints` csv file. This example is based on the PDB structure [1AC5](https://www.rcsb.org/structure/1ac5). The predicted structrue (colored, glycans in purple and orange, protein in green) from this script should look like the following when aligned with the ground truth 1AC5 structure (gray):

![glycan example prediction](./output.png)

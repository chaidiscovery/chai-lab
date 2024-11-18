# Working with restraints

Chai-1 uniquely offers the ability to fold complexes with user-specified "restraints" as inputs. These restraints specify inter-chain contacts at various resolutions that are used to guide Chai-1 in folding the complex. In the following, we describe the file format that Chai-1 expects when working with restraints, as well as an example of how folding might improve with a minimal set of these restraints specified.

## Restraints file format

Restraints to Chai-1 are specified as a table in `.csv` format, where we specify the chain and optionally the residue index for the two partners in a given restraint, along with some information regarding distances, confidence, and metadata. The following example shows what this file might look like:

| restraint_id | chainA | res_idxA | chainB | res_idxB | connection_type | confidence | min_distance_angstrom | max_distance_angstrom | comment |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| restraint0 | A | R84 | C | G7 | contact | 1.0 | 0.0 | 22.0 | toy example |
| restraint1 | C |  | A | S18 | pocket | 1.0 | 0.0 | 11.0 | toy known residue-chain interaction |

The first non-header row `restraint0` shows an example of a `contact` restraint. Chai-1 defines a `contact` restraint to be between two residues in two distinct chains; the distances indicate an upper bound on how far apart we expect those two residues to be. 

The second non-header row `restraint1` shows an example of a `pocket` restraint, which species that a chain is in contact with a specific residue in another chain. Note that compared to `contact` restraints, this is a "coarser" and "assymetric" restraint, as _any_ residue in the first chain can be in contact with the specified token in the second chain. This is notationally indicated by having `res_idxA` being left empty. The aforementioned `contact` restraint, on the other hand, specifies an exact residue for both chains. 

A few other notes:
- In `res_idxA` and `res_idxB` fields, we expect a concatenation of a residue and its 1-based index. For example, if your chain is comprised of the residues `ARNDRA` and you want to specify a contact on the `D` residue, the input should be `D4`. The code internally checks that the given residue at the given index matches the input sequence, and will throw an error if they do not match. This redundancy primarily helps avoid cases where we misspecify positions in constraints.
- Chains `chainA` and `chainB` are assigned identifiers in alphabetical A-Z order following the order that chains are specified in the input. For example, if you wish to specify an interaction between the first and four chains in your input, you'd use the letters `A` and `D`. 
- You may specify a mixture of `contact` and `pocket` restraints.
- The `restraint_id` column must be unique.
- The `comment` field is not read as an input and is included for user convenience.
- The `confidence` and `min_distance_angstrom` fields are _currently_ not used by the model, and are included in this format for future-proofing.

## Example

As an example, consider the PDB structure [7SYZ](https://www.rcsb.org/structure/7SYZ). Without restraints provided, Chai-1 does not accurately predict the interface between the viral protein domain and the heavy/light chains (see below table for interface DockQ scores). We provide $n=2$ randomly selected contact restraints based on the experimental ground truth structure, and see that doing so significantly improves the interface DockQ scores between the viral protein and the antibody chains. Providing $n=2$ pocket restraints has a similarly positive effect, though the effect is smaller due to the lower specificity of pocket restraints (see below table).

| Interface | Chai-1 DockQ w/o restraints | Chai-1 DockQ w/ contact restraints | Chai-1 DockQ w/ pocket restraints |
| --- | --- | --- | --- |
| antibody-light | 0.020 | 0.400 | 0.273 | 
| antibody-heavy | 0.011 | 0.274 | 0.204 | 
| heavy-light | 0.789 | 0.712 | 0.719 |

Code for running this example is given in `examples/restraints/predict_with_restraints.py` along with the restraint files used to produce these results. For this example, DockQ scores are obtained by evaluting the single top-ranked output from the Chai-1 model run without MSAs.
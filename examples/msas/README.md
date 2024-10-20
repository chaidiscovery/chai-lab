# Adding MSA evolutionary information

While Chai-1 performs very well in "single-sequence mode," it can also be given additional evolutionary information to further improve performance. As in other folding methods, this evolutionary information is provided in the form of a multiple sequence alignment (MSA). 

## The `.aligned.pqt` file format

The easiest way to provide MSA information to Chai-1 is through the `.aligned.pqt` file format that we have defined. This file can be thought of as an augmented `a3m` file, and is essentially a dataframe saved in parquet format with the following four (required) columns:

- `sequence` - contains alignment hits in the `a3m` sequence format. 
- `source_database` - contains information regarding the database that the sequence came from. This can be set to one of `uniprot`, `uniref90`, `bfd_uniclust`, `mgnify`, or `query`. This source is featurized as an input to Chai-1 (and is therefore not just for book keeping!). The `query` key should only occur once in the table as the first row and indicates the query sequence used to construct the hits. If your alignments come from a database not included in these four options, it's probably a good idea to experiment with setting a source to get best results; `uniref90` should be a good "catchall" choice in general. 
- `pairing_key` - a string that indicates a "key" for how alignments for different sequences in a complex should be "paired" up with each other (similarly to AlphaFold-Multimer and AlphaFold 3) to capture evolutionary information across chains. Pairing is done across all sequences with the same pairing key. Therefore, while this pairing key is typically set to a species identifier, you might want to consider setting it to some other key that provides a more or less general grouping by which to match up MSA alignments across different sequences. 
- `comment` - the wild west. You may put whatever string you choose here and it will be ignored; this field is provided primarily for human readability and book keeping. 

This file format offers several advantages over standard `.a3m` files:
- Alignments that come from different databases can be contained in a single file without losing track of which alignment came from where. Note, however, that all rows in this file are assumed to be associated with the same (indicated) query sequence. 
- Allows flexible specification of how sequences across different chains ought to be paired up.

See the following for a toy example of what this table might look like:

| sequence | source_database | pairing_key  | comment                       |
| -------- | --------------- | ------------ | ----------------------------- |
| RKDSS... | query           |              | query sequence                |
| RKDES... | uniref90        |              | A fun sequence from uniref90  |
| RKSES... | uniprot         | Mus musculus | A mouse sequence from uniprot |
| ...      |

We additionally provide code to parse `a3m` files into this format; see `merge_multi_a3m_to_aligned_dataframe` in `chai_lab/data/parsing/msas/aligned_pqt.py`.

### TLDR

Chai-1 uses `.aligned.pqt` files to specify MSAs. These are similar to `a3m` with added columns for source database and pairing key to pair MSAs across different chains. Each `.aligned.pqt` file contains all MSAs for a single query sequence. 

## From `.aligned.pqt` to `MSAContext`

By default, the `run_inference` example inference code we provide assumes that all MSAs required for a prediction are stored in a specified folder. Each file corresponds to the all alignments for a given sequence, and filenames are specified by the hash of their sequence (this filename is inferred using code in `chai_lab/data/parsing/msas/aligned_pqt.py`). During inference, the script tries to find `<HASH>.aligned.pqt` files in that folder (one file per unique chain sequence) and loads in a `MSAContext` for each MSA it can find; see `chai_lab/data/dataset/msas/load.py` for details. 

## Putting it all together

To demonstrate how these pieces tie together, we provide `aligned.pqt` files containing MSAs for the example in `examples/predict_structure.py` under the `examples/msas` folder. Inference can be run using these example MSAs by providing the path to this folder as an additional argument to `run_inference` as follows:

```python
candidates = run_inference(
    ...
    msa_directory=Path("examples/msas"),
    ...
)
```

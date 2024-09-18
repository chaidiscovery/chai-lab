from pathlib import Path

import numpy as np
import torch

from chai_lab.chai1 import run_inference

# We use fasta-like format for inputs.
# - each entity encodes protein, ligand, RNA or DNA.
# - each entity is labeled with unique name;
# - ligands are encoded with SMILES; modified residues encoded like AAA(SEP)AAA

# Modify example below to your problem.

example_fasta = """
>protein|example-of-long-protein
AGSHSMRYFSTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRGEPRAPWVEQEGPEYWDRETQKYKRQAQTDRVSLRNLRGYYNQSEAGSHTLQWMFGCDLGPDGRLLRGYDQSAYDGKDYIALNEDLRSWTAADTAAQITQRKWEAAREAEQRRAYLEGTCVEWLRRYLENGKETLQRAEHPKTHVTHHPVSDHEATLRCWALGFYPAEITLTWQWDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGEEQRYTCHVQHEGLPEPLTLRWEP
>protein|example-of-short-protein
AIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM
>protein|example-of-peptide
GAAL
>ligand|and-example-for-ligand-encoded-as-smiles
CCCCCCCCCCCCCC(=O)O
""".strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

output_dir = Path("/tmp/outputs")

candidates = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    # 'default' setup
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device=torch.device("cuda:0"),
    use_esm_embeddings=True,
)

cif_paths = candidates.cif_paths
scores = [rd.aggregate_score for rd in candidates.ranking_data]


# Load pTM, ipTM, pLDDTs and clash scores for sample 2
scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))

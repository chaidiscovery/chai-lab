import datetime
from pathlib import Path

import numpy as np
import torch

import chai_lab
from chai_lab.chai1 import run_inference

# We use fasta-like format for inputs.
# - each entity encodes protein, ligand, RNA or DNA
# - each entity is labeled with unique name;
# - ligands are encoded with SMILES; modified residues encoded like AAA(SEP)AAA

# Example given below, just modify it


example_fasta = (
    """
>protein|name=example-of-long-protein
"""
    + "A" * 1700
).strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

output_dir = (
    Path(__file__)
    .parent.parent.joinpath(
        f"outputs/{chai_lab.__version__}/{datetime.datetime.now().isoformat()}"
    )
    .absolute()
)


if False:
    from chai_lab.dispatch import MemSourceTrackingDispatchMode

    with MemSourceTrackingDispatchMode():
        pass
if True:
    x = torch.zeros(85174583296 // 2, dtype=torch.uint8, device=torch.device("cuda:0"))
    # x = torch.zeros(
    #     85174583296 // 2 - 1_000_000_000,
    #     dtype=torch.uint8,
    #     device=torch.device("cuda:0"),
    # )
    # assert 0 == 1
    candidates = run_inference(
        fasta_file=fasta_path,
        output_dir=output_dir,
        # 'default' setup
        num_trunk_recycles=2,
        num_diffn_timesteps=2,
        seed=42,
        device="cuda:0",
        use_esm_embeddings=True,
        low_memory=True,
    )

cif_paths = candidates.cif_paths
scores = [rd.aggregate_score for rd in candidates.ranking_data]


# Load pTM, ipTM, pLDDTs and clash scores for sample 2
scores = np.load(output_dir.joinpath("scores.model_idx_2.npz"))

import logging
from pathlib import Path

import numpy as np
import torch

from chai_lab.chai1 import run_inference

logging.basicConfig(level=logging.INFO)

example_fasta = """
>protein|7WSL_1-light
DIQLTQSPSFLSAYVGDRVTITCKASQDVGTAVAWYQQKPGKAPKLLIYWASTLHTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQHYSSYPWTFGQGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC
>protein|7WSL_2-PDL1
DRPWNPPTFSPALLVVTEGDNATFTCSFSNTSESFVLNWYRMSPSNQTDKLAAFPEDRSQPGQDSRFRVTQLPNGRDFHMSVVRARRNDSGTYLCGAISLAPKAQIKESLRAELRVTERRAELEHHHHHH
>protein|7WSL_3-heavy
EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYDMSWVRQAPGKGLEWVSTISGGGSYTYYQDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCASPYYAMDYWGQGTTVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHHHHHH
""".strip()

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(example_fasta)

output_dir = Path("/tmp/outputs")

candidates = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    constraint_path="/workspaces/chai-lab/examples/constraints/7WSL.contact.constraints",
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

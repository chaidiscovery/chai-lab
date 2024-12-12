import logging
import shutil
from pathlib import Path

from chai_lab.chai1 import run_inference

logging.basicConfig(level=logging.INFO)

# Inference expects an empty directory; enforce this
output_dir = Path("/workspaces/chai-lab/tmp/outputs")
if output_dir.exists():
    logging.warning(f"Removing old output directory: {output_dir}")
    shutil.rmtree(output_dir)
output_dir.mkdir(exist_ok=True, parents=True)

candidates = run_inference(
    fasta_file=Path(__file__).with_name("1ac5.fasta"),
    output_dir=output_dir,
    constraint_path=Path(__file__).with_name("1ac5.restraints"),
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device="cuda:0",
    use_esm_embeddings=True,
)

cif_paths = candidates.cif_paths

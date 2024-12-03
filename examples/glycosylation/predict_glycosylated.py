import logging
import shutil
from pathlib import Path

from chai_lab.chai1 import run_inference

logging.basicConfig(level=logging.INFO)

glycosylated_fasta = """
>protein|1AC5
LPSSEEYKVAYELLPGLSEVPDPSNIPQMHAGHIPLRSEDADEQDSSDLEYFFWKFTNNDSNGNVDRPLIIWLNGGPGCSSMDGALVESGPFRVNSDGKLYLNEGSWISKGDLLFIDQPTGTGFSVEQNKDEGKIDKNKFDEDLEDVTKHFMDFLENYFKIFPEDLTRKIILSGESYAGQYIPFFANAILNHNKFSKIDGDTYDLKALLIGNGWIDPNTQSLSYLPFAMEKKLIDESNPNFKHLTNAHENCQNLINSASTDEAAHFSYQECENILNLLLSYTRESSQKGTADCLNMYNFNLKDSYPSCGMNWPKDISFVSKFFSTPGVIDSLHLDSDKIDHWKECTNSVGTKLSNPISKPSIHLLPGLLESGIEIVLFNGDKDLICNNKGVLDTIDNLKWGGIKGFSDDAVSFDWIHKSKSTDDSEEFSGYVKYDRNLTFVSVYNASHMVPFDKSLVSRGIVDIYSNDVMIIDNNGKNVMITT
>glycan|two-sugar
NAG(1-4 NAG)
>glycan|one-sugar
NAG
"""

fasta_path = Path("/tmp/example.fasta")
fasta_path.write_text(glycosylated_fasta)

# Inference expects an empty directory; enforce this
output_dir = Path("/workspaces/chai-lab/tmp/outputs")
if output_dir.exists():
    logging.warning(f"Removing old output directory: {output_dir}")
    shutil.rmtree(output_dir)
output_dir.mkdir(exist_ok=True, parents=True)

candidates = run_inference(
    fasta_file=fasta_path,
    output_dir=output_dir,
    constraint_path=Path(__file__).with_name("bonds.restraints"),
    # 'default' setup
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device="cuda:0",
    use_esm_embeddings=True,
)

cif_paths = candidates.cif_paths

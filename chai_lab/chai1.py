# %%
import math
from dataclasses import dataclass
from pathlib import Path

import torch
import torch.export
from einops import einsum, rearrange, repeat
from torch import Tensor
from tqdm import tqdm

from chai_lab.data.collate.collate import Collate
from chai_lab.data.collate.utils import AVAILABLE_MODEL_SIZES
from chai_lab.data.dataset.all_atom_feature_context import (
    MAX_MSA_DEPTH,
    MAX_NUM_TEMPLATES,
    AllAtomFeatureContext,
)
from chai_lab.data.dataset.constraints.constraint_context import ConstraintContext
from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.embeddings.esm import get_esm_embedding_context
from chai_lab.data.dataset.inference_dataset import load_chains_from_raw, read_inputs
from chai_lab.data.dataset.msas.msa_context import MSAContext
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.templates.context import TemplateContext
from chai_lab.data.features.feature_factory import FeatureFactory
from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.atom_element import AtomElementOneHot
from chai_lab.data.features.generators.atom_name import AtomNameOneHot
from chai_lab.data.features.generators.base import EncodingType
from chai_lab.data.features.generators.blocked_atom_pair_distances import (
    BlockedAtomPairDistances,
    BlockedAtomPairDistogram,
)
from chai_lab.data.features.generators.docking import DockingConstraintGenerator
from chai_lab.data.features.generators.esm_generator import ESMEmbeddings
from chai_lab.data.features.generators.identity import Identity
from chai_lab.data.features.generators.is_cropped_chain import ChainIsCropped
from chai_lab.data.features.generators.missing_chain_contact import MissingChainContact
from chai_lab.data.features.generators.msa import (
    IsPairedMSAGenerator,
    MSADataSourceGenerator,
    MSADeletionMeanGenerator,
    MSADeletionValueGenerator,
    MSAFeatureGenerator,
    MSAHasDeletionGenerator,
    MSAProfileGenerator,
)
from chai_lab.data.features.generators.ref_pos import RefPos
from chai_lab.data.features.generators.relative_chain import RelativeChain
from chai_lab.data.features.generators.relative_entity import RelativeEntity
from chai_lab.data.features.generators.relative_sep import RelativeSequenceSeparation
from chai_lab.data.features.generators.relative_token import RelativeTokenSeparation
from chai_lab.data.features.generators.residue_type import ResidueType
from chai_lab.data.features.generators.structure_metadata import (
    IsDistillation,
    TokenBFactor,
    TokenPLDDT,
)
from chai_lab.data.features.generators.templates import (
    TemplateDistogramGenerator,
    TemplateMaskGenerator,
    TemplateResTypeGenerator,
    TemplateUnitVectorGenerator,
)
from chai_lab.data.features.generators.token_dist_restraint import (
    TokenDistanceRestraint,
)
from chai_lab.data.features.generators.token_pair_pocket_restraint import (
    TokenPairPocketRestraint,
)
from chai_lab.data.io.pdb_utils import write_pdbs_from_outputs
from chai_lab.model.diffusion_schedules import InferenceNoiseSchedule
from chai_lab.model.utils import center_random_augmentation
from chai_lab.ranking.frames import get_frames_and_mask
from chai_lab.ranking.rank import SampleRanking, rank
from chai_lab.utils.paths import chai1_component
from chai_lab.utils.tensor_utils import move_data_to_device, set_seed, und_self
from chai_lab.utils.typing import Float, typecheck


class UnsupportedInputError(RuntimeError):
    pass


def load_exported(comp_key: str, device: torch.device) -> torch.nn.Module:
    local_path = chai1_component(comp_key)
    exported_program = torch.export.load(local_path)
    return exported_program.module().to(device)


# %%
# Create feature factory

feature_generators = dict(
    RelativeSequenceSeparation=RelativeSequenceSeparation(sep_bins=None),
    RelativeTokenSeparation=RelativeTokenSeparation(r_max=32),
    RelativeEntity=RelativeEntity(),
    RelativeChain=RelativeChain(),
    ResidueType=ResidueType(
        min_corrupt_prob=0.0,
        max_corrupt_prob=0.0,
        num_res_ty=32,
        key="token_residue_type",
    ),
    ESMEmbeddings=ESMEmbeddings(),  # TODO: this can probably be the identity
    BlockedAtomPairDistogram=BlockedAtomPairDistogram(),
    InverseSquaredBlockedAtomPairDistances=BlockedAtomPairDistances(
        transform="inverse_squared",
        encoding_ty=EncodingType.IDENTITY,
    ),
    AtomRefPos=RefPos(),
    AtomRefCharge=Identity(
        key="inputs/atom_ref_charge",
        ty=FeatureType.ATOM,
        dim=1,
        can_mask=False,
    ),
    AtomRefMask=Identity(
        key="inputs/atom_ref_mask",
        ty=FeatureType.ATOM,
        dim=1,
        can_mask=False,
    ),
    AtomRefElement=AtomElementOneHot(max_atomic_num=128),
    AtomNameOneHot=AtomNameOneHot(),
    TemplateMask=TemplateMaskGenerator(),
    TemplateUnitVector=TemplateUnitVectorGenerator(),
    TemplateResType=TemplateResTypeGenerator(),
    TemplateDistogram=TemplateDistogramGenerator(),
    TokenDistanceRestraint=TokenDistanceRestraint(
        include_probability=0.0,
        size=0.33,
        min_dist=6.0,
        max_dist=30.0,
        num_rbf_radii=6,
    ),
    DockingConstraintGenerator=DockingConstraintGenerator(
        include_probability=0.0,
        structure_dropout_prob=0.75,
        chain_dropout_prob=0.75,
    ),
    TokenPairPocketRestraint=TokenPairPocketRestraint(
        size=0.33,
        include_probability=0.0,
        min_dist=6.0,
        max_dist=20.0,
        coord_noise=0.0,
        num_rbf_radii=6,
    ),
    MSAProfile=MSAProfileGenerator(),
    MSADeletionMean=MSADeletionMeanGenerator(),
    IsDistillation=IsDistillation(),
    TokenBFactor=TokenBFactor(include_prob=0.0),
    TokenPLDDT=TokenPLDDT(include_prob=0.0),
    ChainIsCropped=ChainIsCropped(),
    MissingChainContact=MissingChainContact(contact_threshold=6.0),
    MSAOneHot=MSAFeatureGenerator(),
    MSAHasDeletion=MSAHasDeletionGenerator(),
    MSADeletionValue=MSADeletionValueGenerator(),
    IsPairedMSA=IsPairedMSAGenerator(),
    MSADataSource=MSADataSourceGenerator(),
)
feature_factory = FeatureFactory(feature_generators)

# %%
# Config


class DiffusionConfig:
    S_churn: float = 80
    S_tmin: float = 4e-4
    S_tmax: float = 80.0
    S_noise: float = 1.003
    sigma_data: float = 16.0
    second_order: bool = True


# %%
# Input validation


def raise_if_too_many_tokens(n_actual_tokens: int):
    if n_actual_tokens > max(AVAILABLE_MODEL_SIZES):
        raise UnsupportedInputError(
            f"Too many tokens in input: {n_actual_tokens} > {max(AVAILABLE_MODEL_SIZES)}. "
            "Please limit the length of the input sequence."
        )


def raise_if_too_many_templates(n_actual_templates: int):
    if n_actual_templates > MAX_NUM_TEMPLATES:
        raise UnsupportedInputError(
            f"Too many templates in input: {n_actual_templates} > {MAX_NUM_TEMPLATES}. "
            "Please limit the number of templates."
        )


def raise_if_msa_too_deep(msa_depth: int):
    if msa_depth > MAX_MSA_DEPTH:
        raise UnsupportedInputError(
            f"MSA to deep: {msa_depth} > {MAX_MSA_DEPTH}. "
            "Please limit the MSA depth."
        )


# %%
# Inference logic


@torch.no_grad()
def run_inference(
    fasta_file: Path,
    output_dir: Path,
    use_esm_embeddings: bool = True,
    # expose some params for easy tweaking
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 2,
    seed: int | None = None,
    device: torch.device | None = None,
) -> list[Path]:
    # Prepare inputs
    assert fasta_file.exists(), fasta_file
    fasta_inputs = read_inputs(fasta_file, length_limit=None)
    assert len(fasta_inputs) > 0, "No inputs found in fasta file"

    # Load structure context
    chains = load_chains_from_raw(fasta_inputs)
    contexts = [c.structure_context for c in chains]
    merged_context = AllAtomStructureContext.merge(contexts)
    n_actual_tokens = merged_context.num_tokens
    raise_if_too_many_tokens(n_actual_tokens)

    # Load MSAs
    msa_context = MSAContext.create_empty(
        n_tokens=n_actual_tokens,
        depth=MAX_MSA_DEPTH,
    )
    main_msa_context = MSAContext.create_empty(
        n_tokens=n_actual_tokens,
        depth=MAX_MSA_DEPTH,
    )

    # Load templates
    template_context = TemplateContext.empty(
        n_tokens=n_actual_tokens,
        n_templates=MAX_NUM_TEMPLATES,
    )

    # Load ESM embeddings
    if use_esm_embeddings:
        embedding_context = get_esm_embedding_context(chains, device=device)
    else:
        embedding_context = EmbeddingContext.empty(n_tokens=n_actual_tokens)

    # Constraints
    constraint_context = ConstraintContext.empty()

    # Build final feature context
    feature_context = AllAtomFeatureContext(
        chains=chains,
        structure_context=merged_context,
        msa_context=msa_context,
        main_msa_context=main_msa_context,
        template_context=template_context,
        embedding_context=embedding_context,
        constraint_context=constraint_context,
    )

    output_paths, scores, ranking_data = run_folding_on_context(
        feature_context,
        output_dir=output_dir,
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
        device=device,
    )
    return output_paths


def _bin_centers(min_bin: float, max_bin: float, no_bins: int) -> Tensor:
    return torch.linspace(min_bin, max_bin, 2 * no_bins + 1)[1::2]


@typecheck
@dataclass(frozen=True)
class ConfidenceScores:
    # Predicted aligned error(PAE)
    pae: Float[Tensor, "bs num_tokens num_tokens"]

    # Predicted distance error (PDE)
    pde: Float[Tensor, "bs num_tokens num_tokens"]

    # Predicted local distance difference test (pLDDT)
    plddt: Float[Tensor, "bs num_tokens"]


@torch.no_grad()
def run_folding_on_context(
    feature_context: AllAtomFeatureContext,
    output_dir: Path,
    # expose some params for easy tweaking
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    seed: int | None = None,
    device: torch.device | None = None,
) -> tuple[list[Path], ConfidenceScores, list[SampleRanking]]:
    """
    Function for in-depth explorations.
    User completely controls folding inputs.
    """
    # Set seed
    if seed is not None:
        set_seed([seed])

    if device is None:
        device = torch.device("cuda:0")

    ##
    ## Validate inputs
    ##

    n_actual_tokens = feature_context.structure_context.num_tokens
    raise_if_too_many_tokens(n_actual_tokens)
    raise_if_too_many_templates(feature_context.template_context.num_templates)
    raise_if_msa_too_deep(feature_context.msa_context.depth)
    raise_if_msa_too_deep(feature_context.main_msa_context.depth)

    ##
    ## Prepare batch
    ##

    # Collate inputs into batch
    collator = Collate(
        feature_factory=feature_factory,
        num_key_atoms=128,
        num_query_atoms=32,
    )

    feature_contexts = [feature_context]
    batch_size = len(feature_contexts)
    batch = collator(feature_contexts)
    batch = move_data_to_device(batch, device=device)

    # Get features and inputs from batch
    features = {name: feature for name, feature in batch["features"].items()}
    inputs = batch["inputs"]
    block_indices_h = inputs["block_atom_pair_q_idces"]
    block_indices_w = inputs["block_atom_pair_kv_idces"]
    atom_single_mask = inputs["atom_exists_mask"]
    atom_token_indices = inputs["atom_token_index"].long()
    token_single_mask = inputs["token_exists_mask"]
    token_pair_mask = und_self(token_single_mask, "b i, b j -> b i j")
    token_reference_atom_index = inputs["token_ref_atom_index"]
    atom_within_token_index = inputs["atom_within_token_index"]
    msa_mask = inputs["msa_mask"]
    template_input_masks = und_self(
        inputs["template_mask"], "b t n1, b t n2 -> b t n1 n2"
    )
    block_atom_pair_mask = inputs["block_atom_pair_mask"]

    ##
    ## Load exported models
    ##

    # Model is size-specific
    model_size = min(x for x in AVAILABLE_MODEL_SIZES if n_actual_tokens <= x)

    feature_embedding = load_exported(f"{model_size}/feature_embedding.pt2", device)
    token_input_embedder = load_exported(
        f"{model_size}/token_input_embedder.pt2", device
    )
    trunk = load_exported(f"{model_size}/trunk.pt2", device)
    diffusion_module = load_exported(f"{model_size}/diffusion_module.pt2", device)
    confidence_head = load_exported(f"{model_size}/confidence_head.pt2", device)

    ##
    ## Run the features through the feature embedder
    ##

    embedded_features = feature_embedding.forward(**features)
    token_single_input_feats = embedded_features["TOKEN"]
    token_pair_input_feats, token_pair_structure_input_features = embedded_features[
        "TOKEN_PAIR"
    ].chunk(2, dim=-1)
    atom_single_input_feats, atom_single_structure_input_features = embedded_features[
        "ATOM"
    ].chunk(2, dim=-1)
    block_atom_pair_input_feats, block_atom_pair_structure_input_feats = (
        embedded_features["ATOM_PAIR"].chunk(2, dim=-1)
    )
    template_input_feats = embedded_features["TEMPLATES"]
    msa_input_feats = embedded_features["MSA"]

    ##
    ## Run the inputs through the token input embedder
    ##

    token_input_embedder_outputs: tuple[Tensor, ...] = token_input_embedder.forward(
        token_single_input_feats=token_single_input_feats,
        token_pair_input_feats=token_pair_input_feats,
        atom_single_input_feats=atom_single_input_feats,
        block_atom_pair_feat=block_atom_pair_input_feats,
        block_atom_pair_mask=block_atom_pair_mask,
        block_indices_h=block_indices_h,
        block_indices_w=block_indices_w,
        atom_single_mask=atom_single_mask,
        atom_token_indices=atom_token_indices,
    )
    token_single_initial_repr, token_single_structure_input, token_pair_initial_repr = (
        token_input_embedder_outputs
    )

    ##
    ## Run the input representations through the trunk
    ##

    # Recycle the representations by feeding the output back into the trunk as input for
    # the subsequent recycle
    token_single_trunk_repr = token_single_initial_repr
    token_pair_trunk_repr = token_pair_initial_repr
    for _ in tqdm(range(num_trunk_recycles), desc="Trunk recycles"):
        (token_single_trunk_repr, token_pair_trunk_repr) = trunk.forward(
            token_single_trunk_initial_repr=token_single_initial_repr,
            token_pair_trunk_initial_repr=token_pair_initial_repr,
            token_single_trunk_repr=token_single_trunk_repr,  # recycled
            token_pair_trunk_repr=token_pair_trunk_repr,  # recycled
            msa_input_feats=msa_input_feats,
            msa_mask=msa_mask,
            template_input_feats=template_input_feats,
            template_input_masks=template_input_masks,
            token_single_mask=token_single_mask,
            token_pair_mask=token_pair_mask,
        )

    ##
    ## Denoise the trunk representation by passing it through the diffusion module
    ##

    def _denoise(atom_pos: Tensor, sigma: Tensor, s: int) -> Tensor:
        atom_noised_coords = rearrange(
            atom_pos, "(b s) ... -> b s ...", s=s
        ).contiguous()
        noise_sigma = repeat(sigma, " -> b s", b=batch_size, s=s)
        return diffusion_module.forward(
            token_single_initial_repr=token_single_structure_input.float(),
            token_pair_initial_repr=token_pair_structure_input_features.float(),
            token_single_trunk_repr=token_single_trunk_repr.float(),
            token_pair_trunk_repr=token_pair_trunk_repr.float(),
            atom_single_input_feats=atom_single_structure_input_features.float(),
            atom_block_pair_input_feats=block_atom_pair_structure_input_feats.float(),
            atom_single_mask=atom_single_mask,
            atom_block_pair_mask=block_atom_pair_mask,
            token_single_mask=token_single_mask,
            block_indices_h=block_indices_h,
            block_indices_w=block_indices_w,
            atom_noised_coords=atom_noised_coords.float(),
            noise_sigma=noise_sigma.float(),
            atom_token_indices=atom_token_indices,
        )

    num_diffn_samples = 5  # Fixed at export time
    inference_noise_schedule = InferenceNoiseSchedule(
        s_max=DiffusionConfig.S_tmax,
        s_min=4e-4,
        p=7.0,
        sigma_data=DiffusionConfig.sigma_data,
    )
    sigmas = inference_noise_schedule.get_schedule(
        device=device, num_timesteps=num_diffn_timesteps
    )
    gammas = torch.where(
        (sigmas >= DiffusionConfig.S_tmin) & (sigmas <= DiffusionConfig.S_tmax),
        min(DiffusionConfig.S_churn / num_diffn_timesteps, math.sqrt(2) - 1),
        0.0,
    )

    sigmas_and_gammas = list(zip(sigmas[:-1], sigmas[1:], gammas[:-1]))

    # Initial atom positions
    _, num_atoms = atom_single_mask.shape
    atom_pos = sigmas[0] * torch.randn(
        batch_size * num_diffn_samples, num_atoms, 3, device=device
    )

    for sigma_curr, sigma_next, gamma_curr in tqdm(
        sigmas_and_gammas, desc="Diffusion steps"
    ):
        # Center coords
        atom_pos = center_random_augmentation(
            atom_pos,
            atom_single_mask=repeat(
                atom_single_mask,
                "b a -> (b s) a",
                s=num_diffn_samples,
            ),
        )

        # Alg 2. lines 4-6
        noise = DiffusionConfig.S_noise * torch.randn(
            atom_pos.shape, device=atom_pos.device
        )
        sigma_hat = sigma_curr + gamma_curr * sigma_curr
        atom_pos_noise = (sigma_hat**2 - sigma_curr**2).clamp_min(1e-6).sqrt()
        atom_pos_hat = atom_pos + noise * atom_pos_noise

        # Lines 7-8
        denoised_pos = _denoise(
            atom_pos=atom_pos_hat,
            sigma=sigma_hat,
            s=num_diffn_samples,
        )
        d_i = (atom_pos_hat - denoised_pos) / sigma_hat
        atom_pos = atom_pos_hat + (sigma_next - sigma_hat) * d_i

        # Lines 9-11
        if sigma_next != 0 and DiffusionConfig.second_order:  # second order update
            denoised_pos = _denoise(
                atom_pos,
                sigma=sigma_next,
                s=num_diffn_samples,
            )
            d_i_prime = (atom_pos - denoised_pos) / sigma_next
            atom_pos = atom_pos + (sigma_next - sigma_hat) * ((d_i_prime + d_i) / 2)

    ##
    ## Run the confidence model
    ##

    confidence_outputs: list[tuple[Tensor, ...]] = [
        confidence_head.forward(
            token_single_input_repr=token_single_initial_repr,
            token_single_trunk_repr=token_single_trunk_repr,
            token_pair_trunk_repr=token_pair_trunk_repr,
            token_single_mask=token_single_mask,
            atom_single_mask=atom_single_mask,
            atom_coords=atom_pos[s : s + 1],
            token_reference_atom_index=token_reference_atom_index,
            atom_token_index=atom_token_indices,
            atom_within_token_index=atom_within_token_index,
        )
        for s in range(num_diffn_samples)
    ]

    pae_logits = torch.cat(
        [x[0] for x in confidence_outputs],
    )
    pde_logits = torch.cat(
        [x[1] for x in confidence_outputs],
    )
    plddt_logits = torch.cat(
        [x[2] for x in confidence_outputs],
    )

    assert atom_pos.shape[0] == num_diffn_samples
    assert pae_logits.shape[0] == num_diffn_samples

    def softmax_einsum_and_cpu(
        logits: Tensor, bin_mean: Tensor, pattern: str
    ) -> Tensor:
        # utility to compute score from bin logits
        res = einsum(
            logits.float().softmax(dim=-1), bin_mean.to(logits.device), pattern
        )
        return res.to(device="cpu")

    token_mask_1d = rearrange(token_single_mask, "1 b -> b")

    pae_scores = softmax_einsum_and_cpu(
        pae_logits[:, token_mask_1d, :, :][:, :, token_mask_1d, :],
        _bin_centers(0.0, 32.0, 64),
        "b n1 n2 d, d -> b n1 n2",
    )

    pde_scores = softmax_einsum_and_cpu(
        pde_logits[:, token_mask_1d, :, :][:, :, token_mask_1d, :],
        _bin_centers(0.0, 32.0, 64),
        "b n1 n2 d, d -> b n1 n2",
    )

    plddt_scores_atom = softmax_einsum_and_cpu(
        plddt_logits,
        _bin_centers(0, 1, plddt_logits.shape[-1]),
        "b a d, d -> b a",
    )

    # converting to per-token
    [mask] = atom_single_mask.cpu()
    [indices] = atom_token_indices.cpu()

    def avg_1d(x):
        n = torch.bincount(indices[mask], weights=x[mask])
        d = torch.bincount(indices[mask]).clamp(min=1)
        return n / d

    plddt_scores = torch.stack([avg_1d(x) for x in plddt_scores_atom])

    confidence_scores = ConfidenceScores(
        pae=pae_scores,
        pde=pde_scores,
        plddt=plddt_scores,
    )

    ##
    ## Write the outputs
    ##

    output_paths: list[Path] = []
    ranking_data: list[SampleRanking] = []

    for idx in range(num_diffn_samples):
        # trunk_sample_idx = preds["trunk_sample_index"][idx].item()
        trunk_sample_idx = 0
        out_basename = f"pred.model_trunk_{trunk_sample_idx}_idx_{idx}.pdb"
        pdb_out_path = output_dir / out_basename

        print(f"Writing output to {pdb_out_path}")

        # use 0-100 scale for pLDDT in pdb outputs
        scaled_plddt_scores_per_atom = 100 * plddt_scores_atom[idx : idx + 1]

        write_pdbs_from_outputs(
            coords=atom_pos[idx : idx + 1],
            bfactors=scaled_plddt_scores_per_atom,
            output_batch=move_data_to_device(inputs, torch.device("cpu")),
            write_path=pdb_out_path,
        )
        output_paths.append(pdb_out_path)

        _, valid_frames_mask = get_frames_and_mask(
            atom_pos[idx : idx + 1],
            inputs["token_asym_id"],
            inputs["token_residue_index"],
            inputs["token_backbone_frame_mask"],
            inputs["token_centre_atom_index"],
            inputs["token_exists_mask"],
            inputs["atom_exists_mask"],
            inputs["token_backbone_frame_index"],
            inputs["atom_token_index"],
        )

        ranking_data.append(
            rank(
                atom_pos[idx : idx + 1],
                atom_mask=inputs["atom_exists_mask"],
                atom_token_index=inputs["atom_token_index"],
                token_exists_mask=inputs["token_exists_mask"],
                token_asym_id=inputs["token_asym_id"],
                token_entity_type=inputs["token_entity_type"],
                token_valid_frames_mask=valid_frames_mask,
                lddt_logits=plddt_logits[idx : idx + 1],
                lddt_bin_centers=_bin_centers(0, 1, plddt_logits.shape[-1]).to(
                    plddt_logits.device
                ),
                pae_logits=pae_logits[idx : idx + 1],
                pae_bin_centers=_bin_centers(0.0, 32.0, 64).to(pae_logits.device),
            )
        )

    return output_paths, confidence_scores, ranking_data

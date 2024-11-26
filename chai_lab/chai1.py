# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.


import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
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
from chai_lab.data.dataset.constraints.restraint_context import (
    RestraintContext,
    load_manual_restraints_for_chai1,
)
from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.embeddings.esm import get_esm_embedding_context
from chai_lab.data.dataset.inference_dataset import load_chains_from_raw, read_inputs
from chai_lab.data.dataset.msas.load import get_msa_contexts
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
from chai_lab.data.features.generators.docking import DockingRestraintGenerator
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
from chai_lab.data.io.cif_utils import outputs_to_cif
from chai_lab.data.parsing.restraints import parse_pairwise_table
from chai_lab.model.diffusion_schedules import InferenceNoiseSchedule
from chai_lab.model.utils import center_random_augmentation
from chai_lab.ranking.frames import get_frames_and_mask
from chai_lab.ranking.rank import SampleRanking, get_scores, rank
from chai_lab.utils.paths import chai1_component
from chai_lab.utils.plot import plot_msa
from chai_lab.utils.tensor_utils import move_data_to_device, set_seed, und_self
from chai_lab.utils.typing import Float, typecheck


class UnsupportedInputError(RuntimeError):
    pass


class ModuleWrapper:
    def __init__(self, jit_module):
        self.jit_module = jit_module

    def forward(self, crop_size: int, **kw):
        return getattr(self.jit_module, f"forward_{crop_size}")(**kw)


def load_exported(comp_key: str, device: torch.device) -> ModuleWrapper:
    torch.jit.set_fusion_strategy([("STATIC", 0), ("DYNAMIC", 0)])
    local_path = chai1_component(comp_key)
    assert isinstance(device, torch.device)
    if device != torch.device("cuda:0"):
        # load on cpu first, then move to device
        return ModuleWrapper(torch.jit.load(local_path, map_location="cpu").to(device))
    else:
        # skip loading on CPU.
        return ModuleWrapper(torch.jit.load(local_path).to(device))


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
        include_probability=1.0,
        size=0.33,
        min_dist=6.0,
        max_dist=30.0,
        num_rbf_radii=6,
    ),
    DockingConstraintGenerator=DockingRestraintGenerator(
        include_probability=0.0,
        structure_dropout_prob=0.75,
        chain_dropout_prob=0.75,
    ),
    TokenPairPocketRestraint=TokenPairPocketRestraint(
        size=0.33,
        include_probability=1.0,
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
            f"MSA too deep: {msa_depth} > {MAX_MSA_DEPTH}. "
            "Please limit the MSA depth."
        )


# %%
# Inference logic


@typecheck
@dataclass(frozen=True)
class StructureCandidates:
    # We provide candidates generated by a model,
    #   with confidence predictions and ranking scores.
    # Predicted structure is a candidate with the highest score.

    # locations of CIF files, one file per candidate
    cif_paths: list[Path]
    # scores for each of candidates + info that was used for scoring.
    ranking_data: list[SampleRanking]
    # iff MSA search was performed, we also save a plot as PDF
    msa_coverage_plot_path: Path | None

    # Predicted aligned error(PAE)
    pae: Float[Tensor, "candidate num_tokens num_tokens"]
    # Predicted distance error (PDE)
    pde: Float[Tensor, "candidate num_tokens num_tokens"]
    # Predicted local distance difference test (pLDDT)
    plddt: Float[Tensor, "candidate num_tokens"]

    def __post_init__(self):
        assert len(self.cif_paths) == len(self.ranking_data)
        assert len(self.cif_paths) == len(self.pae)


@torch.no_grad()
def run_inference(
    fasta_file: Path,
    *,
    output_dir: Path,
    use_esm_embeddings: bool = True,
    msa_directory: Path | None = None,
    constraint_path: Path | str | None = None,
    # expose some params for easy tweaking
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    seed: int | None = None,
    device: torch.device | None = None,
) -> StructureCandidates:
    # Prepare inputs
    assert fasta_file.exists(), fasta_file
    fasta_inputs = read_inputs(fasta_file, length_limit=None)

    assert len(fasta_inputs) > 0, "No inputs found in fasta file"

    for name, count in Counter([inp.entity_name for inp in fasta_inputs]).items():
        if count > 1:
            raise UnsupportedInputError(
                f"{name=} used more than once in inputs. Each entity must have a unique name"
            )

    # Load structure context
    chains = load_chains_from_raw(fasta_inputs)
    merged_context = AllAtomStructureContext.merge(
        [c.structure_context for c in chains]
    )
    n_actual_tokens = merged_context.num_tokens
    raise_if_too_many_tokens(n_actual_tokens)

    # Load MSAs
    if msa_directory is not None:
        msa_context, msa_profile_context = get_msa_contexts(
            chains, msa_directory=msa_directory
        )
    else:
        msa_context = MSAContext.create_empty(
            n_tokens=n_actual_tokens, depth=MAX_MSA_DEPTH
        )
        msa_profile_context = MSAContext.create_empty(
            n_tokens=n_actual_tokens, depth=MAX_MSA_DEPTH
        )
    assert (
        msa_context.num_tokens == merged_context.num_tokens
    ), f"Discrepant tokens in input and MSA: {merged_context.num_tokens} != {msa_context.num_tokens}"

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
    if constraint_path is not None:
        restraint_context = load_manual_restraints_for_chai1(
            chains,
            crop_idces=None,
            provided_constraints=parse_pairwise_table(constraint_path),
        )
    else:
        restraint_context = RestraintContext.empty()

    # Build final feature context
    feature_context = AllAtomFeatureContext(
        chains=chains,
        structure_context=merged_context,
        msa_context=msa_context,
        profile_msa_context=msa_profile_context,
        template_context=template_context,
        embedding_context=embedding_context,
        restraint_context=restraint_context,
    )

    return run_folding_on_context(
        feature_context,
        output_dir=output_dir,
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
        device=device,
    )


def _bin_centers(min_bin: float, max_bin: float, no_bins: int) -> Tensor:
    return torch.linspace(min_bin, max_bin, 2 * no_bins + 1)[1::2]


@torch.no_grad()
def run_folding_on_context(
    feature_context: AllAtomFeatureContext,
    *,
    output_dir: Path,
    # expose some params for easy tweaking
    num_trunk_recycles: int = 3,
    num_diffn_timesteps: int = 200,
    seed: int | None = None,
    device: torch.device | None = None,
) -> StructureCandidates:
    """
    Function for in-depth explorations.
    User completely controls folding inputs.
    """
    # Set seed
    if seed is not None:
        set_seed([seed])

    if device is None:
        device = torch.device("cuda:0")

    # Clear memory
    torch.cuda.empty_cache()

    ##
    ## Validate inputs
    ##

    n_actual_tokens = feature_context.structure_context.num_tokens
    raise_if_too_many_tokens(n_actual_tokens)
    raise_if_too_many_templates(feature_context.template_context.num_templates)
    raise_if_msa_too_deep(feature_context.msa_context.depth)
    # NOTE profile MSA used only for statistics; no depth check

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

    _, _, model_size = msa_mask.shape
    assert model_size in AVAILABLE_MODEL_SIZES

    feature_embedding = load_exported("feature_embedding.pt", device)
    token_input_embedder = load_exported("token_embedder.pt", device)
    trunk = load_exported("trunk.pt", device)
    diffusion_module = load_exported("diffusion_module.pt", device)
    confidence_head = load_exported("confidence_head.pt", device)

    ##
    ## Run the features through the feature embedder
    ##

    embedded_features = feature_embedding.forward(crop_size=model_size, **features)
    token_single_input_feats = embedded_features["TOKEN"]
    token_pair_input_feats, token_pair_structure_input_feats = embedded_features[
        "TOKEN_PAIR"
    ].chunk(2, dim=-1)
    atom_single_input_feats, atom_single_structure_input_feats = embedded_features[
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
        crop_size=model_size,
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
            crop_size=model_size,
        )
    # We won't be using the trunk anymore; remove it from memory
    del trunk
    torch.cuda.empty_cache()

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
            token_pair_initial_repr=token_pair_structure_input_feats.float(),
            token_single_trunk_repr=token_single_trunk_repr.float(),
            token_pair_trunk_repr=token_pair_trunk_repr.float(),
            atom_single_input_feats=atom_single_structure_input_feats.float(),
            atom_block_pair_input_feats=block_atom_pair_structure_input_feats.float(),
            atom_single_mask=atom_single_mask,
            atom_block_pair_mask=block_atom_pair_mask,
            token_single_mask=token_single_mask,
            block_indices_h=block_indices_h,
            block_indices_w=block_indices_w,
            atom_noised_coords=atom_noised_coords.float(),
            noise_sigma=noise_sigma.float(),
            atom_token_indices=atom_token_indices,
            crop_size=model_size,
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

    # We won't be running diffusion anymore
    del diffusion_module
    torch.cuda.empty_cache()

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
            crop_size=model_size,
        )
        for s in range(num_diffn_samples)
    ]

    pae_logits, pde_logits, plddt_logits = [
        torch.cat(single_sample, dim=0)
        for single_sample in zip(*confidence_outputs, strict=True)
    ]

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

    # converting per-atom plddt to per-token
    [mask] = atom_single_mask.cpu()
    [indices] = atom_token_indices.cpu()

    def avg_per_token_1d(x):
        n = torch.bincount(indices[mask], weights=x[mask])
        d = torch.bincount(indices[mask]).clamp(min=1)
        return n / d

    plddt_scores = torch.stack([avg_per_token_1d(x) for x in plddt_scores_atom])

    ##
    ## Write the outputs
    ##
    # Move data to the CPU so we don't hit GPU memory limits
    inputs = move_data_to_device(inputs, torch.device("cpu"))
    atom_pos = atom_pos.cpu()
    plddt_logits = plddt_logits.cpu()
    pae_logits = pae_logits.cpu()

    # Plot coverage of tokens by MSA, save plot
    output_dir.mkdir(parents=True, exist_ok=True)

    if feature_context.msa_context.mask.any():
        msa_plot_path = plot_msa(
            input_tokens=feature_context.structure_context.token_residue_type,
            msa_tokens=feature_context.msa_context.tokens,
            out_fname=output_dir / "msa_depth.pdf",
        )
    else:
        msa_plot_path = None

    cif_paths: list[Path] = []
    ranking_data: list[SampleRanking] = []

    for idx in range(num_diffn_samples):
        ##
        ## Compute ranking scores
        ##

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

        ranking_outputs = rank(
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

        ranking_data.append(ranking_outputs)

        ##
        ## Write output files
        ##

        cif_out_path = output_dir.joinpath(f"pred.model_idx_{idx}.cif")
        aggregate_score = ranking_outputs.aggregate_score.item()
        print(f"Score={aggregate_score:.3f}, writing output to {cif_out_path}   ")

        # use 0-100 scale for pLDDT in pdb outputs
        scaled_plddt_scores_per_atom = 100 * plddt_scores_atom[idx : idx + 1]

        outputs_to_cif(
            coords=atom_pos[idx : idx + 1],
            bfactors=scaled_plddt_scores_per_atom,
            output_batch=inputs,
            write_path=cif_out_path,
            entity_names={
                c.entity_data.entity_id: c.entity_data.entity_name
                for c in feature_context.chains
            },
        )
        cif_paths.append(cif_out_path)

        scores_out_path = output_dir.joinpath(f"scores.model_idx_{idx}.npz")

        np.savez(scores_out_path, **get_scores(ranking_outputs))

    return StructureCandidates(
        cif_paths=cif_paths,
        ranking_data=ranking_data,
        msa_coverage_plot_path=msa_plot_path,
        pae=pae_scores,
        pde=pde_scores,
        plddt=plddt_scores,
    )

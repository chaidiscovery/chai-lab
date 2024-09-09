import logging
from dataclasses import asdict, dataclass
from functools import cached_property, partial

import torch
from torch import Tensor

from chai_lab.utils.tensor_utils import (
    batch_tensorcode_to_string,
    tensorcode_to_string,
)
from chai_lab.utils.typing import Bool, Float, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class AllAtomStructureContext:
    # token-level
    token_residue_type: Int[Tensor, "n_tokens"]
    token_residue_index: Int[Tensor, "n_tokens"]
    token_index: Int[Tensor, "n_tokens"]
    token_centre_atom_index: Int[Tensor, "n_tokens"]
    token_ref_atom_index: Int[Tensor, "n_tokens"]
    token_exists_mask: Bool[Tensor, "n_tokens"]
    token_backbone_frame_mask: Bool[Tensor, "n_tokens"]
    token_backbone_frame_index: Int[Tensor, "n_tokens 3"]
    token_asym_id: Int[Tensor, "n_tokens"]
    token_entity_id: Int[Tensor, "n_tokens"]
    token_sym_id: Int[Tensor, "n_tokens"]
    token_entity_type: Int[Tensor, "n_tokens"]
    token_residue_name: UInt8[Tensor, "n_tokens 8"]
    token_b_factor_or_plddt: Float[Tensor, "n_tokens"]
    # atom-level
    atom_token_index: Int[Tensor, "n_atoms"]
    atom_within_token_index: Int[Tensor, "n_atoms"]  # consistent atom ordering
    atom_ref_pos: Float[Tensor, "n_atoms 3"]
    atom_ref_mask: Bool[Tensor, "n_atoms"]
    atom_ref_element: Int[Tensor, "n_atoms"]
    atom_ref_charge: Int[Tensor, "n_atoms"]
    atom_ref_name: list[str]
    atom_ref_name_chars: Int[Tensor, "n_atoms 4"]
    atom_ref_space_uid: Int[Tensor, "n_atoms"]
    atom_is_not_padding_mask: Bool[Tensor, "n_atoms"]
    # supervision only
    atom_gt_coords: Float[Tensor, "n_atoms 3"]
    atom_exists_mask: Bool[Tensor, "n_atoms"]
    # structure-level
    pdb_id: UInt8[Tensor, "n_tokens 32"]
    # source_pdb_chain_id corresponds to auth_asym_id in pdb
    # can be the same for two different asym_id values
    # (we split protein and ligand for example)
    source_pdb_chain_id: UInt8[Tensor, "n_tokens 4"]
    # subchain_id is label_asym_id in pdb
    # it is assigned by the PDB and separates different
    # chemical entities (protein, ligand)
    # should be a 1-1 mapping to asym_id
    subchain_id: UInt8[Tensor, "n_tokens 4"]
    resolution: Float[Tensor, "1"]
    is_distillation: Bool[Tensor, "1"]
    # symmetric atom swap indices
    symmetries: Int[Tensor, "n_atoms n_symmetries"]

    def __post_init__(self):
        # Resolved residues filter should eliminate PDBs with missing residues, but that
        # we can still have atom_exists mask set to False at every position if we have a
        # bad crop so we log examples with no valid coordinates
        if self.num_atoms > 0 and not torch.any(self.atom_exists_mask):
            pdb_id = tensorcode_to_string(self.pdb_id[0])
            logger.error(f"No valid coordinates found in any atoms for {pdb_id}")

        # Check that atom and token masks are compatible. Anywhere that the atom mask is
        # true, the token mask should also be true
        if self.num_atoms > 0 and not torch.all(
            self.token_exists_mask[self.atom_token_index][self.atom_exists_mask]
        ):
            pdb_id = tensorcode_to_string(self.pdb_id[0])
            logger.error(f"Incompatible masks for {pdb_id}")

    @cached_property
    def residue_names(self) -> list[str]:
        return batch_tensorcode_to_string(self.token_residue_name)

    def pad(
        self,
        n_tokens: int,
        n_atoms: int,
    ) -> "AllAtomStructureContext":
        assert n_tokens >= self.num_tokens
        pad_tokens_func = partial(_pad_func, pad_size=n_tokens - self.num_tokens)

        assert n_atoms >= self.num_atoms
        pad_atoms_func = partial(_pad_func, pad_size=n_atoms - self.num_atoms)

        return AllAtomStructureContext(
            # token-level
            token_residue_type=pad_tokens_func(self.token_residue_type),
            token_residue_index=pad_tokens_func(self.token_residue_index),
            token_index=pad_tokens_func(self.token_index),
            token_centre_atom_index=pad_tokens_func(self.token_centre_atom_index),
            token_ref_atom_index=pad_tokens_func(self.token_ref_atom_index),
            token_exists_mask=pad_tokens_func(self.token_exists_mask),
            token_backbone_frame_mask=pad_tokens_func(self.token_backbone_frame_mask),
            token_backbone_frame_index=torch.cat(
                [
                    pad_tokens_func(self.token_backbone_frame_index[..., i]).unsqueeze(
                        -1
                    )
                    for i in range(3)
                ],
                dim=-1,
            ),
            token_asym_id=pad_tokens_func(self.token_asym_id),
            token_entity_id=pad_tokens_func(self.token_entity_id),
            token_sym_id=pad_tokens_func(self.token_sym_id),
            token_entity_type=pad_tokens_func(self.token_entity_type),
            token_residue_name=pad_tokens_func(self.token_residue_name),
            token_b_factor_or_plddt=pad_tokens_func(self.token_b_factor_or_plddt),
            # atom-level
            atom_token_index=pad_atoms_func(self.atom_token_index),
            atom_within_token_index=pad_atoms_func(self.atom_within_token_index),
            atom_ref_pos=pad_atoms_func(self.atom_ref_pos),
            atom_ref_mask=pad_atoms_func(self.atom_ref_mask),
            atom_ref_element=pad_atoms_func(self.atom_ref_element),
            atom_ref_charge=pad_atoms_func(self.atom_ref_charge),
            atom_ref_name=self.atom_ref_name,
            atom_ref_name_chars=pad_atoms_func(self.atom_ref_name_chars),
            atom_ref_space_uid=pad_atoms_func(self.atom_ref_space_uid, pad_value=-1),
            atom_is_not_padding_mask=pad_atoms_func(self.atom_is_not_padding_mask),
            # supervision-only
            atom_gt_coords=pad_atoms_func(self.atom_gt_coords),
            atom_exists_mask=pad_atoms_func(self.atom_exists_mask),
            # structure-level
            pdb_id=pad_tokens_func(self.pdb_id),
            source_pdb_chain_id=pad_tokens_func(self.source_pdb_chain_id),
            subchain_id=pad_tokens_func(self.subchain_id),
            resolution=self.resolution,
            is_distillation=self.is_distillation,
            symmetries=pad_atoms_func(self.symmetries, pad_value=-1),
        )

    @typecheck
    @classmethod
    def merge(
        cls,
        contexts: list["AllAtomStructureContext"],
    ) -> "AllAtomStructureContext":
        # indexes:
        token_offsets = _exclusive_cum_lengths([x.token_residue_type for x in contexts])
        atom_offsets = _exclusive_cum_lengths([x.atom_token_index for x in contexts])

        atom_token_index = torch.cat(
            [x.atom_token_index + count for x, count in zip(contexts, token_offsets)]
        )

        token_centre_atom_index = torch.cat(
            [
                x.token_centre_atom_index + count
                for x, count in zip(contexts, atom_offsets)
            ]
        )
        token_ref_atom_index = torch.cat(
            [x.token_ref_atom_index + count for x, count in zip(contexts, atom_offsets)]
        )
        token_backbone_frame_index = torch.cat(
            [
                x.token_backbone_frame_index + count
                for x, count in zip(contexts, token_offsets)
            ]
        )

        n_tokens = sum(x.num_tokens for x in contexts)
        token_index = torch.arange(n_tokens, dtype=torch.int)

        # re-index the reference space from 0..n_tokens-1.
        zero_indexed_ref_uids = [
            torch.unique_consecutive(x.atom_ref_space_uid, return_inverse=True)[1]
            for x in contexts
        ]

        ref_space_uids_offsets = _exclusive_cum_lengths(
            [x.atom_ref_space_uid for x in contexts]
        )
        atom_ref_space_uid = torch.cat(
            [
                x + count
                for x, count in zip(zero_indexed_ref_uids, ref_space_uids_offsets)
            ],
        )

        # pad symmetric permutations to have same length
        max_symms = max(x.symmetries.shape[-1] for x in contexts)
        padded_symms = [
            torch.nn.functional.pad(
                x.symmetries, (0, max_symms - x.symmetries.shape[-1]), value=-1
            )
            for x in contexts
        ]
        # offset symmetries by number of atoms in each chain
        symm_mask = torch.cat([x >= 0 for x in padded_symms])
        symmetries = torch.cat(padded_symms)
        symmetries = symmetries.masked_fill(~symm_mask, -1)

        return cls(
            # token-level
            token_residue_type=torch.cat([x.token_residue_type for x in contexts]),
            token_residue_index=torch.cat([x.token_residue_index for x in contexts]),
            token_index=token_index,
            token_centre_atom_index=token_centre_atom_index,
            token_ref_atom_index=token_ref_atom_index,
            token_exists_mask=torch.cat([x.token_exists_mask for x in contexts]),
            token_backbone_frame_mask=torch.cat(
                [x.token_backbone_frame_mask for x in contexts]
            ),
            token_backbone_frame_index=token_backbone_frame_index,
            token_asym_id=torch.cat([x.token_asym_id for x in contexts]),
            token_entity_id=torch.cat([x.token_entity_id for x in contexts]),
            token_sym_id=torch.cat([x.token_sym_id for x in contexts]),
            token_entity_type=torch.cat([x.token_entity_type for x in contexts]),
            token_residue_name=torch.cat([x.token_residue_name for x in contexts]),
            token_b_factor_or_plddt=torch.cat(
                [x.token_b_factor_or_plddt for x in contexts]
            ),
            # atom-level
            atom_token_index=atom_token_index,
            atom_within_token_index=torch.cat(
                [x.atom_within_token_index for x in contexts]
            ),
            atom_ref_pos=torch.cat([x.atom_ref_pos for x in contexts]),
            atom_ref_mask=torch.cat([x.atom_ref_mask for x in contexts]),
            atom_ref_element=torch.cat([x.atom_ref_element for x in contexts]),
            atom_ref_charge=torch.cat([x.atom_ref_charge for x in contexts]),
            atom_ref_name=[x for context in contexts for x in context.atom_ref_name],
            atom_ref_name_chars=torch.cat([x.atom_ref_name_chars for x in contexts]),
            atom_ref_space_uid=atom_ref_space_uid,
            atom_is_not_padding_mask=torch.cat(
                [x.atom_is_not_padding_mask for x in contexts]
            ),
            # supervision only
            atom_gt_coords=torch.cat([x.atom_gt_coords for x in contexts]),
            atom_exists_mask=torch.cat([x.atom_exists_mask for x in contexts]),
            # structure-level
            pdb_id=torch.cat([x.pdb_id for x in contexts]),
            source_pdb_chain_id=torch.cat([x.source_pdb_chain_id for x in contexts]),
            subchain_id=torch.cat([x.subchain_id for x in contexts]),
            resolution=torch.max(
                torch.stack([x.resolution for x in contexts]), 0
            ).values,
            is_distillation=torch.max(
                torch.stack([x.is_distillation for x in contexts]), 0
            ).values,
            symmetries=symmetries,
        )

    def to(self, device: torch.device | str) -> "AllAtomStructureContext":
        dict_ = {
            k: v.to(device) if torch.is_tensor(v) else v
            for k, v in asdict(self).items()
        }
        return AllAtomStructureContext(**dict_)

    @property
    def num_tokens(self) -> int:
        (n_tokens,) = self.token_index.shape
        return n_tokens

    @property
    def num_atoms(self) -> int:
        (n_atoms,) = self.atom_token_index.shape
        return n_atoms

    def to_dict(self) -> dict[str, torch.Tensor]:
        return asdict(self)


def _pad_func(x: Tensor, pad_size: int, pad_value: float | None = None) -> Tensor:
    sizes = [0, 0] * (x.ndim - 1) + [0, pad_size]
    return torch.nn.functional.pad(x, sizes, value=pad_value)


def _exclusive_cum_lengths(tensors: list[Int[Tensor, "n"]]):
    lengths = torch.tensor([t.shape[0] for t in tensors])
    cum_lengths = torch.cumsum(lengths, dim=0).roll(1, 0)
    cum_lengths[0] = 0
    return cum_lengths

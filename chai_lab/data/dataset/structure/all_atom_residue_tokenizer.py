# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import dataclass
from itertools import chain

import torch
from einops import repeat
from torch import Tensor

from chai_lab.data.dataset.structure import utils
from chai_lab.data.dataset.structure.all_atom_structure_context import (
    AllAtomStructureContext,
)
from chai_lab.data.dataset.structure.bond_utils import (
    get_atom_covalent_bond_pairs_from_glycan_string,
)
from chai_lab.data.dataset.structure.utils import (
    backbone_atoms_all_present,
    backbone_atoms_indices,
    get_centre_atom_index,
    get_reference_atom_index,
)
from chai_lab.data.parsing.structure.all_atom_entity_data import AllAtomEntityData
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.data.parsing.structure.residue import ConformerData, Residue
from chai_lab.data.residue_constants import standard_residue_pdb_codes
from chai_lab.data.sources.rdkit import (
    RefConformerGenerator,
    conformer_data_to_rdkit_mol,
)
from chai_lab.utils.tensor_utils import string_to_tensorcode, unique_indexes
from chai_lab.utils.typing import Bool, Float, Int, typecheck

logger = logging.getLogger(__name__)


# jaxtyping on residue-level objects is extremely slow.
@dataclass(frozen=True)
class TokenSpan:
    restype: Int[Tensor, "n_tokens"]
    residue_index: Int[Tensor, "n_tokens"]
    centre_atom_index: Int[Tensor, "n_tokens"]
    reference_atom_index: Int[Tensor, "n_tokens"]
    backbone_frame_mask: Bool[Tensor, "n_tokens"]
    backbone_frame_index: Int[Tensor, "n_tokens 3"]
    atom_gt_coords: Float[Tensor, "n_atoms 3"]
    atom_exists_mask: Bool[Tensor, "n_atoms"]
    atom_token_index: Int[Tensor, "n_atoms"]
    ref_pos: Float[Tensor, "n_atoms 3"]
    ref_mask: Bool[Tensor, "n_atoms"]
    ref_element: Int[Tensor, "n_atoms"]
    ref_charge: Int[Tensor, "n_atoms"]
    atom_names: list[str]
    # Consistent atom ordering witin each token
    atom_within_token_indices: Int[Tensor, "n_atoms"]
    residue_names: list[str]
    symmetries: Int[Tensor, "n_atoms n_symm"]
    b_factor_or_plddt: Float[Tensor, "n_tokens"]

    @classmethod
    def concatenate(cls, spans: list["TokenSpan"]) -> "TokenSpan":
        # offset bond indices:
        tokens_per_span = torch.tensor([span.restype.shape[0] for span in spans])
        token_count = torch.cumsum(tokens_per_span, dim=0).roll(1, 0)
        token_count[0] = 0

        # offsets indices of centre atoms:
        atoms_per_span = torch.tensor(
            [span.atom_exists_mask.shape[0] for span in spans]
        )
        atom_offsets = torch.cumsum(atoms_per_span, dim=0).roll(1, 0)
        atom_offsets[0] = 0

        centre_atom_index = torch.cat(
            [
                span.centre_atom_index + offset
                for span, offset in zip(spans, atom_offsets)
            ]
        )
        reference_atom_index = torch.cat(
            [
                span.reference_atom_index + offset
                for span, offset in zip(spans, atom_offsets)
            ]
        )

        atom_token_index = (
            torch.cumsum(
                torch.cat([x.atom_token_index for x in spans]),
                dim=0,
                dtype=torch.int,
            )
            - 1
        )
        backbone_frame_index = torch.cat(
            [
                span.backbone_frame_index + offset
                for span, offset in zip(spans, atom_offsets)
            ]
        )

        # concatenate symmetric permutations at the atom level
        # make sure that trailing shape is the same
        # NOTE: we store the *local* permutation indices, not the global ones
        # i.e. the permutation indices are relative to the residue
        atom_symms = [span.symmetries for span in spans]
        max_symms = max(x.shape[-1] for x in atom_symms)
        atom_symms = [
            torch.nn.functional.pad(x, (0, max_symms - x.shape[-1]), value=-1)
            for x in atom_symms
        ]
        return cls(
            restype=torch.cat([x.restype for x in spans]),
            residue_index=torch.cat([x.residue_index for x in spans]),
            centre_atom_index=centre_atom_index,
            reference_atom_index=reference_atom_index,
            backbone_frame_mask=torch.cat([x.backbone_frame_mask for x in spans]),
            backbone_frame_index=backbone_frame_index,
            atom_gt_coords=torch.cat([x.atom_gt_coords for x in spans]),
            atom_exists_mask=torch.cat([x.atom_exists_mask for x in spans]),
            atom_token_index=atom_token_index,
            ref_pos=torch.cat([x.ref_pos for x in spans]),
            ref_mask=torch.cat([x.ref_mask for x in spans]),
            ref_element=torch.cat([x.ref_element for x in spans]),
            ref_charge=torch.cat([x.ref_charge for x in spans]),
            atom_names=list(chain.from_iterable([x.atom_names for x in spans])),
            atom_within_token_indices=torch.cat(
                [x.atom_within_token_indices for x in spans]
            ),
            residue_names=list(chain.from_iterable([x.residue_names for x in spans])),
            symmetries=torch.cat(atom_symms, dim=0),
            b_factor_or_plddt=torch.cat([x.b_factor_or_plddt for x in spans]),
        )


class AllAtomResidueTokenizer:
    ref_conformer_generator: RefConformerGenerator

    def __init__(self, ref_conformer_generator: RefConformerGenerator):
        self.ref_conformer_generator = ref_conformer_generator

    def tokenize_residue(
        self,
        residue: Residue,
        entity_type: EntityType,
    ) -> TokenSpan | None:
        ref_conformer_data = self._get_ref_conformer_data(residue)
        if ref_conformer_data.num_atoms == 0:
            # avoid dealing with empty tensors in downstream processing
            # this should only happen when residue is sole hydrogen
            # or when residue code is not in CCD dictionary and
            # the residue has 0 coords in the PDB structure
            logger.warning(
                f"skipping residue {residue.name} {residue.label_seq} as reference conformer has 0 heavy atoms"
            )
            return None

        # Keep only the atoms from the ground truth conformer that are present in
        # reference conformer.
        #
        # If we don't have a reference conformer, we fall back to using the ground truth
        # conformer names, i.e. we keep all atoms in the ground truth conformer.
        # When a true conformer data is not provided, use reference conformer directly
        gt_conformer_data = residue.conformer_data

        if gt_conformer_data is not None:
            atom_gt_coords, atom_exists_mask = gt_conformer_data.gather_atom_positions(
                ref_conformer_data.atom_names
            )
        else:
            atom_gt_coords = ref_conformer_data.position
            atom_exists_mask = torch.ones(
                atom_gt_coords.shape[0], dtype=torch.bool, device=atom_gt_coords.device
            )

        # Tokenization is by residue if it is a standard amino acid or standard
        # nucleotide; all ligands and all modified residues are tokenized per atom.
        tokenize_fn = (
            self._tokenize_per_residue
            if (
                residue.name in standard_residue_pdb_codes
                and entity_type != EntityType.LIGAND
                and entity_type != EntityType.MANUAL_GLYCAN
            )
            else self._tokenize_per_atom
        )

        return tokenize_fn(
            restype=torch.tensor([residue.restype], dtype=torch.int),
            residue_index=torch.tensor([residue.residue_index], dtype=torch.int),
            atom_gt_coords=atom_gt_coords,
            atom_exists_mask=atom_exists_mask,
            ref_pos=ref_conformer_data.position,
            ref_mask=torch.ones_like(atom_exists_mask),
            ref_element=ref_conformer_data.element,
            ref_charge=ref_conformer_data.charge,
            atom_names=ref_conformer_data.atom_names,
            residue_name=residue.name,
            bonds=ref_conformer_data.bonds,
            symmetries=ref_conformer_data.symmetries,
            b_factor_or_plddt=torch.tensor([residue.b_factor_or_plddt]),
        )

    @staticmethod
    def filter_atom_symmetries(
        symmetries: Int[Tensor, "n_atoms n_symm"],
        atom_exists_mask: Bool[Tensor, "n_atoms"],
    ) -> Int[Tensor, "n_atoms filtered_n_symm"]:
        n_atoms, _ = symmetries.shape

        # Create a mask for non-trivial symmetries
        atom_indices = torch.arange(n_atoms).unsqueeze(-1)
        non_trivial_symmetries = (symmetries >= 0) & (symmetries != atom_indices)

        masked_atoms = ~atom_exists_mask.unsqueeze(-1)

        # Check if any of the masked-out atoms have non-trivial symmetries
        violations = torch.any(masked_atoms & non_trivial_symmetries, dim=1)

        # If any invalid symmetries are found, replace with identity permutation
        if torch.any(violations):
            return atom_indices

        # Otherwise, return the original symmetries
        return symmetries

    # jaxtyping on residue-level objects is very slow,
    # use for debug only
    # @typecheck
    def _tokenize_per_residue(
        self,
        restype: Int[Tensor, "n_tokens"],
        residue_index: Int[Tensor, "n_tokens"],
        atom_gt_coords: Float[Tensor, "n_atoms 3"],
        atom_exists_mask: Bool[Tensor, "n_atoms"],
        ref_pos: Float[Tensor, "n_atoms 3"],
        ref_mask: Bool[Tensor, "n_atoms"],
        ref_element: Int[Tensor, "n_atoms"],
        ref_charge: Int[Tensor, "n_atoms"],
        atom_names: list[str],
        residue_name: str,
        bonds: list[tuple[int, int]],
        symmetries: Int[Tensor, "n_atoms n_symm"],
        b_factor_or_plddt: Float[Tensor, "n_tokens"],
    ) -> TokenSpan:
        centre_atom_index = get_centre_atom_index(
            atom_names,
            residue_name,
        )
        reference_atom_index = get_reference_atom_index(
            atom_names,
            residue_name,
        )
        backbone_frame_mask = backbone_atoms_all_present(
            atom_names,
            residue_name,
        )
        backbone_indices = backbone_atoms_indices(atom_names, residue_name).unsqueeze(0)

        # to 1 token
        atom_token_index = torch.zeros_like(atom_exists_mask, dtype=torch.int)
        atom_token_index[0] = 1

        residue_names = [residue_name]

        # Find atom ordering; these should always be available because per residue
        # tokenization works only on standard residues.
        atom_within_token_index = atom_names_to_atom37_indices(
            atom_names=atom_names,
            residue_name=residue_name,
        )

        return TokenSpan(
            restype=restype,
            residue_index=residue_index,
            centre_atom_index=centre_atom_index,
            reference_atom_index=reference_atom_index,
            backbone_frame_mask=backbone_frame_mask,
            backbone_frame_index=backbone_indices,
            atom_gt_coords=atom_gt_coords,
            atom_exists_mask=atom_exists_mask,
            atom_token_index=atom_token_index,
            ref_pos=ref_pos,
            ref_mask=ref_mask,
            ref_element=ref_element,
            ref_charge=ref_charge,
            atom_names=atom_names,
            atom_within_token_indices=atom_within_token_index,
            residue_names=residue_names,
            symmetries=self.filter_atom_symmetries(symmetries, atom_exists_mask),
            b_factor_or_plddt=b_factor_or_plddt,
        )

    # jaxtyping on residue-level objects is very slow,
    # use for debug only
    # @typecheck
    def _tokenize_per_atom(
        self,
        restype: Int[Tensor, "n_tokens"],
        residue_index: Int[Tensor, "n_tokens"],
        atom_gt_coords: Float[Tensor, "n_atoms 3"],
        atom_exists_mask: Bool[Tensor, "n_atoms"],
        ref_pos: Float[Tensor, "n_atoms 3"],
        ref_mask: Bool[Tensor, "n_atoms"],
        ref_element: Int[Tensor, "n_atoms"],
        ref_charge: Int[Tensor, "n_atoms"],
        atom_names: list[str],
        residue_name: str,
        bonds: list[tuple[int, int]],
        symmetries: Int[Tensor, "n_atoms n_symm"],
        b_factor_or_plddt: Float[Tensor, "n_tokens"],
    ) -> TokenSpan:
        # to n_atoms tokens
        n_atoms = atom_gt_coords.shape[0]
        restype = repeat(restype, "1 -> a", a=n_atoms)
        residue_index = repeat(residue_index, "1 -> a", a=n_atoms)
        b_factor_or_plddt = repeat(b_factor_or_plddt, "1 -> a", a=n_atoms)

        # centre of the token is the first and only atom in each token
        # when tokenizing per-atom
        centre_atom_index = torch.arange(n_atoms, dtype=torch.int)
        reference_atom_index = torch.arange(n_atoms, dtype=torch.int)
        backbone_frame_mask = torch.zeros((n_atoms,), dtype=torch.bool)
        backbone_indices = (
            torch.arange(n_atoms, dtype=torch.int).unsqueeze(1).expand(-1, 3)
        )

        atom_token_index = torch.ones_like(atom_exists_mask, dtype=torch.int)

        residue_names = [residue_name] * n_atoms

        # Each atom is alone in its own token
        atom_within_token_index = torch.zeros(n_atoms, dtype=torch.int)

        return TokenSpan(
            restype=restype,
            residue_index=residue_index,
            centre_atom_index=centre_atom_index,
            reference_atom_index=reference_atom_index,
            backbone_frame_mask=backbone_frame_mask,
            backbone_frame_index=backbone_indices,
            atom_gt_coords=atom_gt_coords,
            atom_exists_mask=atom_exists_mask,
            atom_token_index=atom_token_index,
            ref_pos=ref_pos,
            ref_mask=ref_mask,
            ref_element=ref_element,
            ref_charge=ref_charge,
            atom_names=atom_names,
            atom_within_token_indices=atom_within_token_index,
            residue_names=residue_names,
            symmetries=self.filter_atom_symmetries(symmetries, atom_exists_mask),
            b_factor_or_plddt=b_factor_or_plddt,
        )

    def tokenize_entity(
        self, entity_data: AllAtomEntityData
    ) -> AllAtomStructureContext | None:
        return self.tokenize_entities([entity_data])[0]

    def tokenize_entities(
        self,
        entities_data: list[AllAtomEntityData],
    ) -> list[AllAtomStructureContext | None]:
        sym_ids = _make_sym_ids([x.entity_id for x in entities_data])

        return [
            self._tokenize_entity(
                entity_data,
                chain_id=idx + 1,
                sym_id=sym_id,
            )
            for idx, (entity_data, sym_id) in enumerate(zip(entities_data, sym_ids))
        ]

    def _tokenize_entity(
        self,
        entity_data: AllAtomEntityData,
        chain_id: int = 1,
        sym_id: int = 1,
    ) -> AllAtomStructureContext | None:
        tokenized_residues = [
            self.tokenize_residue(residue, entity_data.entity_type)
            for residue in entity_data.residues
        ]

        valid_residues: list[TokenSpan] = [
            x for x in tokenized_residues if x is not None
        ]
        if len(valid_residues) == 0:
            logger.warning(
                f"Got no residues for entity {entity_data.entity_id} with residues {entity_data.residues}"
            )
            return None

        tokens = TokenSpan.concatenate(valid_residues)

        num_tokens = tokens.restype.shape[0]
        token_index = torch.arange(num_tokens, dtype=torch.int)

        # mask indicating if a token has >=1 atom with known coordinates
        token_exists_mask = (tokens.atom_token_index == token_index[..., None]).sum(
            dim=-1
        ) > 0

        # checks on atom mask and positions:
        # max 1 atom per-example has zero coordinates
        if (
            torch.sum(
                torch.all(tokens.atom_gt_coords[tokens.atom_exists_mask] == 0, dim=-1)
            )
            > 1
        ):
            raise ValueError(
                f"Zero coordinates found in unmasked atoms for {entity_data.pdb_id}"
            )

        # construct asym_id, entity_id, sym_id
        asym_id = chain_id
        entity_id = entity_data.entity_id

        # Create unique ids to identify atoms which belong to same residue in same chain
        # here assume we featurize a single chain
        atom_residue_index = torch.gather(
            tokens.residue_index,
            dim=0,
            index=tokens.atom_token_index.long(),
        )

        atom_ref_space_uid = atom_residue_index

        residue_names = tokens.residue_names

        match entity_data.entity_type:
            case EntityType.PROTEIN:
                if tokens.residue_index[0].item() != 0:
                    logger.error(
                        f"Protein residue index should start at zero, {entity_data}"
                    )

                if not torch.all(torch.diff(tokens.residue_index) <= 1):
                    logger.error(
                        f"Protein residue index should be contiguous (no gaps), {entity_data}"
                    )

                _, unique_indices = unique_indexes(tokens.residue_index)
                res_seq = [residue_names[i.item()] for i in unique_indices]
                if res_seq != entity_data.full_sequence:
                    logger.error(
                        f"Protein residue names should match entity data full sequence, {entity_data}"
                    )

        return AllAtomStructureContext(
            # token-level
            token_residue_type=tokens.restype,
            token_residue_index=tokens.residue_index,
            token_centre_atom_index=tokens.centre_atom_index,
            token_ref_atom_index=tokens.reference_atom_index,
            token_index=token_index,
            token_exists_mask=token_exists_mask,
            token_backbone_frame_mask=tokens.backbone_frame_mask,
            token_backbone_frame_index=tokens.backbone_frame_index,
            token_asym_id=_id_to_token_tensor(asym_id, num_tokens),
            token_entity_id=_id_to_token_tensor(entity_id, num_tokens),
            token_sym_id=_id_to_token_tensor(sym_id, num_tokens),
            token_entity_type=entity_type_to_tensor(
                entity_data.entity_type,
                num_tokens,
            ),
            # token res name is padded to 8 characters
            token_residue_name=torch.stack(
                [string_to_tensorcode(x, 8) for x in residue_names],
                dim=0,
            ),
            token_b_factor_or_plddt=tokens.b_factor_or_plddt,
            # atom-level
            atom_token_index=tokens.atom_token_index,
            atom_within_token_index=tokens.atom_within_token_indices,
            atom_ref_pos=tokens.ref_pos,
            atom_ref_mask=tokens.ref_mask,
            atom_ref_element=tokens.ref_element,
            atom_ref_charge=tokens.ref_charge,
            atom_ref_name=tokens.atom_names,
            atom_ref_name_chars=_atom_names_to_tensor(tokens.atom_names),
            atom_ref_space_uid=atom_ref_space_uid,
            atom_is_not_padding_mask=torch.ones_like(
                tokens.atom_exists_mask,
                dtype=torch.bool,
            ),
            # supervision only
            atom_gt_coords=tokens.atom_gt_coords,
            atom_exists_mask=tokens.atom_exists_mask,
            # structure-only
            pdb_id=repeat(
                # PDB ids are only 4 characters long, but AFDB ids can be longer
                string_to_tensorcode(entity_data.pdb_id, pad_to_length=32),
                "length -> num_tokens length",
                num_tokens=num_tokens,
            ),
            source_pdb_chain_id=repeat(
                string_to_tensorcode(entity_data.source_pdb_chain_id, pad_to_length=4),
                "length -> num_tokens length",
                num_tokens=num_tokens,
            ),
            subchain_id=repeat(
                string_to_tensorcode(entity_data.subchain_id, pad_to_length=4),
                "length -> num_tokens length",
                num_tokens=num_tokens,
            ),
            resolution=torch.tensor(
                [entity_data.resolution],
                dtype=torch.float32,
            ),
            is_distillation=torch.tensor(
                [entity_data.is_distillation],
                dtype=torch.bool,
            ),
            symmetries=tokens.symmetries,
            atom_covalent_bond_indices=get_atom_covalent_bond_pairs_from_glycan_string(
                glycan_string=(
                    entity_data.original_record
                    if entity_data.entity_type == EntityType.MANUAL_GLYCAN
                    else ""
                ),
                token_residue_index=tokens.residue_index,
                atom_ref_name=tokens.atom_names,
            ),
        )

    def _get_ref_conformer_data(self, residue: Residue) -> ConformerData:
        """
        Returns the reference conformer data for the residue. We determine the reference
        conformer according to the following logic:
        1. conformer_generator is available and a reference
           conformer exists for the residue name => we return the cached reference
           conformer via the conformer generator
        2. conformer_generator is available and a smiles is given for the residue =>
            we generate a reference conformer using Rdkit via the conformer generator
        3. conformer_generator is available, the reference conformer can't be
            found and no smiles is given => we convert the Residue to an RDKit molecule
            and load full conformer data with the residue's atom positions as coordinates.
        4. conformer generator is not available => we set reference conformer to
            the ground truth conformer
        """
        # The reference conformer tells us:
        # - which atoms we should expect in this ligand / residue, and how many of them
        # - what are the ideal coordinates of these atoms if the ligand or residue was
        #   assembled alone in the void
        ref_conformer = self.ref_conformer_generator.get(residue.name)

        if ref_conformer is not None:
            if residue.name in standard_residue_pdb_codes:
                return ref_conformer
            else:
                return ref_conformer.center_random_augment()

        # When we can't find a reference conformer, and a smiles is given,
        # generate a reference conformer using rdkit
        if residue.smiles is not None:
            logger.info(
                f"Generating ref conformer for {residue.name}, {residue.smiles}"
            )
            return self.ref_conformer_generator.generate(residue.smiles)

        # When we can't find a reference conformer, attempt to use the ground
        # truth conformer data as the reference conformer.
        logger.warning(
            f"No reference conformer found for residue {residue.name},"
            "using training example conformer"
        )
        assert residue.conformer_data is not None

        try:
            # Rather than just setting the reference conformer to the ground truth, we
            # make a fake RDKit molecule from the ground truth data and then convert
            # back into a conformer data so that we can extract inter-atom aymmetries
            # bond and info
            rdkit_mol = conformer_data_to_rdkit_mol(residue.conformer_data)
            gt_conformer = RefConformerGenerator._load_ref_conformer_from_rdkit(
                rdkit_mol
            )
        except Exception as e:
            # Occasionally _load_ref_conformer_from_rdkit fails on unknown ligands e.g.
            # rdkit.Chem.rdchem.AtomValenceException's can be raised or ValueError:
            # can't infer bonds for Ligand. due to inexact connectivity.
            logger.warning(
                f"Caught error for {residue.name=} while loading reference conformer "
                f"from RDKit, {(type(e).__name__)}. Using ground truth conformer instead."
            )
            gt_conformer = residue.conformer_data

        return gt_conformer.center_random_augment()


@typecheck
def _atom_names_to_tensor(atom_names: list[str]) -> Int[Tensor, "n_atoms 4"]:
    ords = torch.tensor(
        [[ord(c) - 32 for c in atom_name.ljust(4, " ")] for atom_name in atom_names],
        dtype=torch.int,
    )
    return ords[:, :4]


@typecheck
def _id_to_token_tensor(id: int, num_tokens: int) -> Int[Tensor, "n"]:
    return id * torch.ones((num_tokens,), dtype=torch.int)


@typecheck
def entity_type_to_tensor(entity_type: EntityType, num_tokens: int) -> Int[Tensor, "n"]:
    return torch.full((num_tokens,), fill_value=entity_type.value, dtype=torch.int)


def _make_sym_ids(entity_ids_per_chain: list[int]) -> list[int]:
    entities_dict: dict[int, int] = dict()
    sym_ids = []

    for entity_id in entity_ids_per_chain:
        sym_id = entities_dict.get(entity_id, 0)
        sym_ids.append(sym_id)
        entities_dict[entity_id] = sym_id + 1

    return sym_ids


def atom_names_to_atom37_indices(
    atom_names: list[str], residue_name: str
) -> Int[Tensor, "n_atoms"]:
    """
    Returns a tensor of indices into the token-level atom names.
    """
    # Proteins use the atom37 ordering and indexing
    # nucleotides use the 36 atom ordering and indexing
    # - DNA is written as DA DG DC DT
    # - RNA is given as A G C U

    precomputed_idces = utils.atom_37_atom_indices()

    if residue_name == "UNK":
        retval = torch.arange(len(atom_names), dtype=torch.int)

    elif residue_name in standard_residue_pdb_codes:
        idx = [precomputed_idces[(residue_name, atom_name)] for atom_name in atom_names]
        retval = torch.tensor(idx, dtype=torch.int)
    else:
        raise ValueError(
            f"Unknown residue name {residue_name} (atom names: {atom_names})"
        )

    assert retval.max() <= 36, f"Out of bounds ordering {atom_names} in {residue_name}"
    return retval

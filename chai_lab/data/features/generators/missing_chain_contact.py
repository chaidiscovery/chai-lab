# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import torch
from torch import Tensor

from chai_lab.data.features.feature_type import FeatureType
from chai_lab.data.features.generators.base import EncodingType, FeatureGenerator
from chai_lab.utils.tensor_utils import cdist, und_self
from chai_lab.utils.typing import Bool, Float, Int, typecheck


class MissingChainContact(FeatureGenerator):
    contact_threshold: float

    def __init__(
        self,
        # Use DockQ atom contact cutoff as default
        contact_threshold: float = 6.0,
    ):
        """Token-Level feature that indicates is a chain has no tokens
        in contact with tokens from another chain.
        """
        super().__init__(
            ty=FeatureType.TOKEN,
            can_mask=False,
            encoding_ty=EncodingType.IDENTITY,
            num_classes=1,
            mult=1,
        )
        self.contact_threshold = contact_threshold

    def get_input_kwargs_from_batch(self, batch) -> dict:
        return dict(
            atom_gt_coords=batch["inputs"]["atom_gt_coords"],
            atom_exists_mask=batch["inputs"]["atom_exists_mask"],
            token_exists_mask=batch["inputs"]["token_exists_mask"],
            token_asym_id=batch["inputs"]["token_asym_id"].long(),
            atom_token_index=batch["inputs"]["atom_token_index"].long(),
        )

    @typecheck
    def _generate(
        self,
        atom_gt_coords: Float[Tensor, "b a 3"],
        atom_exists_mask: Bool[Tensor, "b a"],
        token_exists_mask: Bool[Tensor, "b n"],
        token_asym_id: Int[Tensor, "b n"],
        atom_token_index: Int[Tensor, "b a"],
    ) -> Tensor:
        # per-atom asym id
        atom_asym_id = torch.gather(token_asym_id, dim=1, index=atom_token_index.long())
        # compute atom pair distances and mask
        atom_pair_dist = cdist(atom_gt_coords)
        atom_pair_mask = und_self(atom_exists_mask, "b i, b j -> b i j")
        atom_pair_asym_mask = atom_asym_id.unsqueeze(-1) != atom_asym_id.unsqueeze(-2)
        aggregate_mask = (
            atom_pair_mask
            & atom_pair_asym_mask
            & (atom_pair_dist < self.contact_threshold)
        )
        # determine which atoms are in contact with some atom from another chain
        atom_in_contact = aggregate_mask.any(dim=-1)
        # determine if any chain has no atoms in contact with another chain
        chain_contact_features: list[torch.Tensor] = []
        for b in range(atom_gt_coords.shape[0]):
            unique_chain_asyms = torch.unique(token_asym_id[b][token_exists_mask[b]])
            if len(unique_chain_asyms) == 1:
                # monomers are set to have no missing contacts
                chain_contact_features.append(
                    torch.zeros_like(
                        token_asym_id[b].unsqueeze(-1), dtype=torch.float32
                    )
                )
                continue
            unique_asyms_with_contacts = torch.unique(
                atom_asym_id[b][atom_in_contact[b]]
            )
            unique_chain_asyms, unique_asyms_with_contacts = [
                set(x.tolist())
                for x in (unique_chain_asyms, unique_asyms_with_contacts)
            ]
            asyms_without_contacts = torch.tensor(
                list(unique_chain_asyms - unique_asyms_with_contacts)
            )
            # create feature data for this chain
            feat = torch.any(
                token_asym_id[b].unsqueeze(-1) == asyms_without_contacts,
                dim=-1,
                keepdim=True,
            )
            chain_contact_features.append(feat.float())

        # make the feature
        return self.make_feature(torch.stack(chain_contact_features, dim=0))

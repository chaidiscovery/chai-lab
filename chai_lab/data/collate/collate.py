# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import dataclasses
import logging
from typing import Any

import torch

from chai_lab.data.collate.utils import get_pad_sizes
from chai_lab.data.dataset.all_atom_feature_context import AllAtomFeatureContext
from chai_lab.data.features.feature_factory import FeatureFactory
from chai_lab.model.utils import (
    get_block_atom_pair_mask,
    get_qkv_indices_for_blocks,
)
from chai_lab.utils.dict import list_dict_to_dict_list

logger = logging.getLogger(__name__)


@dataclasses.dataclass(frozen=True)
class Collate:
    feature_factory: FeatureFactory
    num_query_atoms: int
    num_key_atoms: int

    def __call__(
        self,
        feature_contexts: list[AllAtomFeatureContext],
    ) -> dict[str, Any]:
        raw_batch = self._collate(feature_contexts)
        prepared_batch = self._post_collate(raw_batch)
        return prepared_batch

    def _collate(
        self,
        feature_contexts: list[AllAtomFeatureContext],
    ) -> dict[str, Any]:
        # Get the pad sizes, finding the max number of tokens/atoms/bonds in the batch.
        pad_sizes = get_pad_sizes([p.structure_context for p in feature_contexts])

        # Pad each feature context to the max sizes
        padded_feature_contexts = [
            feature_context.pad(
                n_tokens=pad_sizes.n_tokens,
                n_atoms=pad_sizes.n_atoms,
            )
            for feature_context in feature_contexts
        ]

        # Convert all the input data into dicts, for each feature context
        inputs_per_context = [e.to_dict() for e in padded_feature_contexts]

        # Stack the dict inputs into a single batch dict, across all feature contexts
        batched_inputs = {
            k: (torch.stack(v, dim=0) if isinstance(v[0], torch.Tensor) else v)
            for k, v in list_dict_to_dict_list(inputs_per_context).items()
        }

        # Make a batch dict
        batch = dict(inputs=batched_inputs)
        return batch

    def _post_collate(self, raw_batch: dict[str, Any]) -> dict[str, Any]:
        """
        takes a list of processed multi-chain systems,
        returns a dictionary with batched tensors to feed in the model forward method
        and any other necessary data for the task/losses
        """
        raw_b_i = raw_batch["inputs"]

        # prepare atom pair block data:
        atom_exists_mask = raw_b_i["atom_exists_mask"]
        block_q_atom_idces, block_kv_atom_idces, kv_mask = get_qkv_indices_for_blocks(
            atom_exists_mask.shape[1],
            self.num_query_atoms,
            self.num_key_atoms,
            atom_exists_mask.device,
        )
        block_atom_pair_mask = get_block_atom_pair_mask(
            atom_single_mask=raw_b_i["atom_ref_mask"],
            q_idx=block_q_atom_idces,
            kv_idx=block_kv_atom_idces,
            kv_is_wrapped_mask=kv_mask,
        )
        raw_b_i |= dict(
            block_atom_pair_q_idces=block_q_atom_idces,
            block_atom_pair_kv_idces=block_kv_atom_idces,
            block_atom_pair_mask=block_atom_pair_mask,
        )

        # Compute features
        raw_batch["features"] = self.feature_factory.generate(raw_batch)

        return raw_batch

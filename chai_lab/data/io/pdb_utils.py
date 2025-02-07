# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from dataclasses import asdict, dataclass
from functools import cached_property

from torch import Tensor

from chai_lab.utils.tensor_utils import tensorcode_to_string
from chai_lab.utils.typing import Bool, Int, UInt8, typecheck

logger = logging.getLogger(__name__)


@typecheck
@dataclass
class PDBContext:
    """Complex (multiple entities) represented as a collection of tensors"""

    token_residue_index: Int[Tensor, "n_tokens"]
    token_asym_id: Int[Tensor, "n_tokens"]
    token_entity_type: Int[Tensor, "n_tokens"]
    token_entity_id: Int[Tensor, "n_tokens"]
    token_residue_names: UInt8[Tensor, "n_tokens 8"]
    token_centre_atom_index: Int[Tensor, "n_tokens"]
    atom_token_index: Int[Tensor, "n_atoms"]
    atom_ref_element: Int[Tensor, "n_atoms"]
    atom_exists_mask: Bool[Tensor, "n_atoms"]
    atom_ref_name_chars: Int[Tensor, "n_atoms 4"]

    @cached_property
    def token_res_names_to_string(self) -> list[str]:
        return [tensorcode_to_string(x) for x in self.token_residue_names.cpu()]

    @cached_property
    def asym_id2entity_type(self) -> dict[int, int]:
        return {
            asym_id: ent_type
            for asym_id, ent_type in zip(
                self.token_asym_id.tolist(), self.token_entity_type.tolist()
            )
            if asym_id != 0
        }


def pdb_context_from_batch(d: dict) -> PDBContext:
    context = PDBContext(
        token_residue_index=d["token_residue_index"][0],
        token_asym_id=d["token_asym_id"][0],
        token_entity_type=d["token_entity_type"][0],
        token_entity_id=d["token_entity_id"][0],
        token_residue_names=d["token_residue_name"][0],
        token_centre_atom_index=d["token_centre_atom_index"][0],
        atom_token_index=d["atom_token_index"][0],
        atom_ref_element=d["atom_ref_element"][0],
        atom_exists_mask=d["atom_exists_mask"][0],
        atom_ref_name_chars=d["atom_ref_name_chars"][0],
    )
    for k, v in asdict(context).items():
        assert v.device.type == "cpu", ("not on cpu:", k, v.device)
    return context

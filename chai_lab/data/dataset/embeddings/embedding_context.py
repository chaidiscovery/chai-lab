# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

from dataclasses import asdict, dataclass

import torch
from torch import Tensor

from chai_lab.utils.typing import Float, typecheck


@typecheck
@dataclass
class EmbeddingContext:
    esm_embeddings: Float[Tensor, "num_tokens d_emb"]

    def __str__(self) -> str:
        return (
            f"{self.__class__.__name__}(esm_embeddings of {self.esm_embeddings.shape})"
        )

    @property
    def num_tokens(self) -> int:
        (num_tokens, _) = self.esm_embeddings.shape
        return num_tokens

    def pad(self, max_tokens: int) -> "EmbeddingContext":
        assert self.num_tokens <= max_tokens

        pad_dims_token = (0, max_tokens - self.num_tokens)
        pad_dims_emb = (0, 0)

        padded_embeddings = torch.nn.functional.pad(
            self.esm_embeddings,
            pad_dims_emb + pad_dims_token,
            value=0,
        )

        return EmbeddingContext(
            esm_embeddings=padded_embeddings,
        )

    def to_dict(self) -> dict[str, torch.Tensor]:
        return asdict(self)

    @classmethod
    def empty(cls, n_tokens: int, d_emb: int = 2560) -> "EmbeddingContext":
        return cls(
            esm_embeddings=torch.zeros(n_tokens, d_emb),
        )

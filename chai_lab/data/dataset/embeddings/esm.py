# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import os
from contextlib import contextmanager

import torch

from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.paths import download_if_not_exists, downloads_path
from chai_lab.utils.typing import typecheck

_esm_model: list = []  # persistent in-process container

os.register_at_fork(after_in_child=lambda: _esm_model.clear())


ESM_URL = "https://chaiassets.com/chai1-inference-depencencies/esm2/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt"


esm_cache_folder = downloads_path.joinpath("esm")


@contextmanager
def esm_model(device):
    """Context transiently keeps ESM model on specified device."""

    local_esm_path = downloads_path.joinpath(
        "esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt"
    )
    download_if_not_exists(ESM_URL, local_esm_path)

    if len(_esm_model) == 0:
        # lazy loading of the model
        if device != torch.device("cuda:0"):
            # load on cpu first, then move to device
            model = torch.jit.load(local_esm_path, map_location="cpu").to(device)
        else:
            # skip loading on CPU.
            model = torch.jit.load(local_esm_path).to(device)

        _esm_model.append(model)

    [model] = _esm_model
    model.to(device)
    model.eval()
    yield model
    model.to("cpu")  # move model back to CPU when done


token_map = {
    "<cls>": 0,
    "<pad>": 1,
    "<eos>": 2,
    "<unk>": 3,
    "L": 4,
    "A": 5,
    "G": 6,
    "V": 7,
    "S": 8,
    "E": 9,
    "R": 10,
    "T": 11,
    "I": 12,
    "D": 13,
    "P": 14,
    "K": 15,
    "Q": 16,
    "N": 17,
    "F": 18,
    "Y": 19,
    "M": 20,
    "H": 21,
    "W": 22,
    "C": 23,
    "X": 24,
    "B": 25,
    "U": 26,
    "Z": 27,
    "O": 28,
    ".": 29,
    "-": 30,
    "<null_1>": 31,
    "<mask>": 32,
}


class DumbTokenizer:
    def __init__(self, token_map: dict[str, int]):
        self.token_map = token_map

    def tokenize(self, text: str) -> list[int]:
        tokens = []
        i = 0
        while i < len(text):
            for token in self.token_map:
                if text.startswith(token, i):
                    tokens.append(self.token_map[token])
                    i += len(token)
                    break
            else:
                raise RuntimeError("Unknown token: " + text[i:])
        return tokens


esm_tokenizer = DumbTokenizer(token_map=token_map)


def _get_esm_contexts_for_sequences(
    prot_sequences: set[str], device
) -> dict[str, EmbeddingContext]:
    if len(prot_sequences) == 0:
        return {}  # skip loading ESM

    seq2embedding_context = {}

    with torch.no_grad():
        with esm_model(device=device) as model:
            for seq in prot_sequences:
                # add bos/eos, tokenize
                token_ids = torch.asarray(esm_tokenizer.tokenize(f"<cls>{seq}<eos>"))
                token_ids = token_ids[None, :].to(device)

                last_hidden_state = model(tokens=token_ids)
                # remove BOS/EOS, back to CPU
                esm_embeddings = last_hidden_state[0, 1:-1].float().to("cpu")
                seq_len, _emb_dim = esm_embeddings.shape
                assert seq_len == len(seq)

                seq2embedding_context[seq] = EmbeddingContext(
                    esm_embeddings=esm_embeddings
                )

    return seq2embedding_context


@typecheck
def get_esm_embedding_context(chains: list[Chain], device) -> EmbeddingContext:
    # device is used for computing, but result is still on CPU

    protein_seq2emb_context = _get_esm_contexts_for_sequences(
        prot_sequences=set(
            chain.entity_data.sequence
            for chain in chains
            if chain.entity_data.entity_type == EntityType.PROTEIN
        ),
        device=device,
    )

    chain_embs = []
    for chain in chains:
        if chain.entity_data.entity_type == EntityType.PROTEIN:
            chain_embs.append(protein_seq2emb_context[chain.entity_data.sequence])
        else:
            # embed non-proteins with zeros
            chain_embs.append(
                EmbeddingContext.empty(n_tokens=chain.structure_context.num_tokens)
            )

    exploded_embs = [
        embedding.esm_embeddings[chain.structure_context.token_residue_index, :]
        for embedding, chain in zip(chain_embs, chains, strict=True)
    ]

    # don't crop any chains during inference
    cropped_embs = exploded_embs

    # if we had to crop, we'd need some logic like below:
    # crop_idces: list[torch.Tensor]
    # cropped_embs = [
    #     embedding[crop_idx, :] for embedding, crop_idx in zip(exploded_embs, crop_idces)
    # ]

    # Merge the embeddings along the tokens dimension (i.e. merge the chains)
    merged_embs = torch.cat(cropped_embs, dim=0)

    return EmbeddingContext(esm_embeddings=merged_embs)

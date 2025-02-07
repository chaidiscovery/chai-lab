# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import os
from contextlib import contextmanager

import torch
from transformers import logging as tr_logging

from chai_lab.data.dataset.embeddings.embedding_context import EmbeddingContext
from chai_lab.data.dataset.structure.chain import Chain
from chai_lab.data.parsing.structure.entity_type import EntityType
from chai_lab.utils.paths import downloads_path
from chai_lab.utils.tensor_utils import move_data_to_device
from chai_lab.utils.typing import typecheck

_esm_model: list = []  # persistent in-process container

os.register_at_fork(after_in_child=lambda: _esm_model.clear())


# unfortunately huggingface complains on pooler layer in ESM being non-initialized.
# Did not find a way to filter specifically that logging message :/
tr_logging.set_verbosity_error()

esm_cache_folder = downloads_path.joinpath("esm")


@contextmanager
def esm_model(model_name: str, device):
    """Context transiently keeps ESM model on specified device."""
    from transformers import EsmModel

    if len(_esm_model) == 0:
        # lazy loading of the model
        _esm_model.append(
            EsmModel.from_pretrained(model_name, cache_dir=esm_cache_folder)
        )

    [model] = _esm_model
    model.to(device)
    model.eval()
    yield model
    model.to("cpu")  # move model back to CPU when done


def _get_esm_contexts_for_sequences(
    prot_sequences: set[str], device
) -> dict[str, EmbeddingContext]:
    if len(prot_sequences) == 0:
        return {}  # skip loading ESM

    # local import, requires huggingface transformers
    from transformers import EsmTokenizer

    model_name = "facebook/esm2_t36_3B_UR50D"
    tokenizer = EsmTokenizer.from_pretrained(model_name, cache_dir=esm_cache_folder)

    seq2embedding_context = {}

    with torch.no_grad():
        with esm_model(model_name=model_name, device=device) as model:
            for seq in prot_sequences:
                inputs = tokenizer(seq, return_tensors="pt")
                inputs = move_data_to_device(dict(**inputs), device=device)
                outputs = model(**inputs)
                # remove BOS/EOS, back to CPU
                esm_embeddings = outputs.last_hidden_state[0, 1:-1].to("cpu")
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

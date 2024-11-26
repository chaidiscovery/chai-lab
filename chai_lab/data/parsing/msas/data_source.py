# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
from enum import Enum

logger = logging.getLogger(__name__)


class MSADataSource(Enum):
    # Special value for the query sequence
    QUERY = "query"

    UNIPROT = "uniprot"
    UNIREF90 = "uniref90"
    BFD = "BFD"
    MGNIFY = "mgnify"
    PAIRED = "paired"
    MAIN = "main"
    BFD_UNICLUST = "bfd_uniclust"
    SINGLETON = "singleton"

    # pad value
    NONE = "none"

    # templates
    PDB70 = "pdb70"

    # ran with 3 jackhmmer iterations (-N=3),
    # higher quality but sloow to generate
    UNIPROT_N3 = "uniprot_n3"
    UNIREF90_N3 = "uniref90_n3"
    MGNIFY_N3 = "mgnify_n3"

    @classmethod
    def get_default_sources(cls):
        return [
            MSADataSource.BFD_UNICLUST,
            MSADataSource.MGNIFY,
            MSADataSource.UNIREF90,
            MSADataSource.UNIPROT,
        ]


def encode_source_to_int(source: MSADataSource) -> int:
    return msa_dataset_source_to_int.get(source, 4)


msa_dataset_source_to_quota: dict[MSADataSource, int] = {
    MSADataSource.UNIREF90: 10_000,
    MSADataSource.UNIPROT: 50_000,
    MSADataSource.BFD_UNICLUST: 1_000_000,  # unlimited
    MSADataSource.BFD: 5000,
    MSADataSource.MGNIFY: 5000,
    MSADataSource.UNIREF90_N3: 10_000,
    MSADataSource.UNIPROT_N3: 50_000,
    MSADataSource.MGNIFY_N3: 5000,
    MSADataSource.PDB70: 5000,
}
msa_dataset_source_to_priority = {
    db: i for i, db in enumerate(msa_dataset_source_to_quota)
}
# query should always go first
msa_dataset_source_to_priority[MSADataSource.QUERY] = -1

# This becomes a feature so changing it might break checkpoint compatibility
msa_dataset_source_to_int = {
    MSADataSource.BFD_UNICLUST: 0,
    MSADataSource.MGNIFY: 1,
    MSADataSource.UNIREF90: 2,
    MSADataSource.UNIPROT: 3,
    MSADataSource.NONE: 4,
    MSADataSource.UNIPROT_N3: 3,
    MSADataSource.UNIREF90_N3: 2,
    MSADataSource.MGNIFY_N3: 1,
    MSADataSource.QUERY: 5,  # in chai-1 remapped to none.
}

database_ids: set[str] = set(x.value for x in MSADataSource)

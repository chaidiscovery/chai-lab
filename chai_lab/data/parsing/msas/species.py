# Copyright (c) 2024 Chai Discovery, Inc.
# This source code is licensed under the Chai Discovery Community License
# Agreement (LICENSE.md) found in the root directory of this source tree.

import hashlib
import logging
import re
from functools import lru_cache

from chai_lab.data.parsing.msas.data_source import MSADataSource

logger = logging.getLogger(__name__)

UNKNOWN_SPECIES = 0


@lru_cache(maxsize=1_000_000)
def stable_hash(str) -> int:
    return int(hashlib.sha256(str.encode("utf-8")).hexdigest()[:7], 16)


def get_tax_names(descriptions: list[str], msa_data_source: MSADataSource) -> list[str]:
    """Extract taxonomy strings from descriptions; empty string if not found."""
    match msa_data_source:
        case MSADataSource.UNIPROT | MSADataSource.UNIPROT_N3:
            # example: ...tein OS=Sphenostylis stenocarpa OX=92480...
            compiled = re.compile(r"OS=([^\s]+ [^\s]+)")
            return [
                match.group(1) if (match := compiled.search(s)) else ""
                for s in descriptions
            ]
        case MSADataSource.UNIREF90:
            compiled = re.compile(r"TaxID=(\d+)")
            return [
                match.group(1) if (match := compiled.search(s)) else ""
                for s in descriptions
            ]
        case _:
            return [""] * len(descriptions)


def get_tax_ids(descriptions: list[str], msa_data_source: MSADataSource) -> list[int]:
    """Extract hash-based tax IDs from descriptions."""
    tax_names = get_tax_names(descriptions, msa_data_source)
    return [
        stable_hash(tax_name) if tax_name else UNKNOWN_SPECIES for tax_name in tax_names
    ]

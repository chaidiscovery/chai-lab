# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import logging
import re

from chai_lab.data.parsing.msas.data_source import MSADataSource

logger = logging.getLogger(__name__)


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

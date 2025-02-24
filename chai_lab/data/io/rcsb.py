# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.
from pathlib import Path

from chai_lab.utils.paths import download_if_not_exists


def download_cif_file(pdb_id: str, directory: Path) -> Path:
    """Download the cif file for the given PDB ID from RCSB into the directory."""
    outfile = directory / f"{pdb_id}.cif.gz"
    source_url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    download_if_not_exists(source_url, outfile)
    assert outfile.exists() and outfile.stat().st_size > 0
    return outfile

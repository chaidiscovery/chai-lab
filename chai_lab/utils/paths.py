# Copyright (c) 2024 Chai Discovery, Inc.
# Licensed under the Apache License, Version 2.0.
# See the LICENSE file for details.

import dataclasses
import logging
import os
import random
from pathlib import Path
from typing import Final

import requests
from filelock import FileLock

# use this path object to specify location
# of anything within repository
repo_root: Final[Path] = Path(__file__).parents[2].absolute()

# weights and helper data is downloaded to CHAI_DOWNLOADS_DIR if provided.
# otherwise we use <repo>/downloads, which is gitignored by default
downloads_path = repo_root.joinpath("downloads")
downloads_path = Path(os.environ.get("CHAI_DOWNLOADS_DIR", downloads_path))


# minimal sanity check in case we start moving things around
assert repo_root.exists()


def download_if_not_exists(http_url: str, path: Path):
    if path.exists():
        return

    with FileLock(path.with_suffix(".download_lock")):
        if path.exists():
            return  # if-lock-if sandwich to download only once
        logging.info(f"downloading {http_url}")
        tmp_path = path.with_suffix(f".download_tmp_{random.randint(10 ** 5, 10**6)}")
        with requests.get(http_url, stream=True) as response:
            response.raise_for_status()  # Check if the request was successful
            # Open a local file with the specified name
            path.parent.mkdir(exist_ok=True, parents=True)
            with tmp_path.open("wb") as file:
                # Download the file in chunks
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:  # Filter out keep-alive new chunks
                        file.write(chunk)
    tmp_path.rename(path)
    assert path.exists()


@dataclasses.dataclass
class Downloadable:
    url: str
    path: Path

    def get_path(self) -> Path:
        # downloads artifact if necessary
        download_if_not_exists(self.url, path=self.path)
        return self.path


cached_conformers = Downloadable(
    url="https://chaiassets.com/chai1-inference-depencencies/conformers_v1.apkl",
    path=downloads_path.joinpath("conformers_v1.apkl"),
)

COMPONENT_URL = (
    "https://chaiassets.com/chai1-inference-depencencies/models_v2/{comp_key}"
)


def chai1_component(comp_key: str) -> Path:
    """
    Downloads exported model, stores in locally in the repo/downloads
    comp_key: e.g. 'trunk.pt'
    """
    assert comp_key.endswith(".pt")
    url = COMPONENT_URL.format(comp_key=comp_key)
    result = downloads_path.joinpath("models_v2", comp_key)
    download_if_not_exists(url, result)

    return result

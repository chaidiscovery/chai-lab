import dataclasses
from pathlib import Path
from typing import Final

import requests

# use this path object to specify location
# of anything within repository
repo_root: Final[Path] = Path(__file__).parents[2].absolute()

# minimal sanity check in case we start moving things around
assert repo_root.exists()


def download(http_url: str, path: Path):
    print(f"downloading {http_url}")
    tmp_path = path.with_suffix(".download_tmp")

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
        if not self.path.exists():
            download(self.url, path=self.path)

        return self.path


cached_conformers = Downloadable(
    url="https://chaiassets.com/chai1-inference-depencencies/conformers.apkl",
    path=repo_root.joinpath("downloads", "conformers.apkl"),
)


def chai1_component(comp_key: str) -> Path:
    """
    Downloads exported model, stores in locally in the repo/downloads
    comp_key: e.g. '384/trunk.pt2'
    """
    assert comp_key.endswith(".pt2")
    url = f"https://chaiassets.com/chai1-inference-depencencies/models/{comp_key}"
    result = repo_root.joinpath("downloads", "models", comp_key)
    if not result.exists():
        download(url, result)

    return result

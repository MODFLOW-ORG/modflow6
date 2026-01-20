import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from warnings import warn

import flopy
import pytest
from conftest import project_root_path
from flaky import flaky
from modflow_devtools.build import meson_build
from modflow_devtools.download import download_and_unzip, get_release
from modflow_devtools.misc import get_ostag

repository = "MODFLOW-ORG/modflow6"
top_bin_path = project_root_path / "bin"


def get_asset_name(asset: dict) -> str:
    ostag = get_ostag()
    name = asset["name"]
    if "win" in ostag:
        return name
    else:
        prefix = name.rpartition("_")[0]
        prefix += f"_{ostag}"
        return f"{prefix}.zip"


@pytest.fixture
def rebuilt_bin_path() -> Path:
    return top_bin_path / "rebuilt"


@pytest.fixture
def downloaded_bin_path() -> Path:
    return top_bin_path / "downloaded"


@flaky(max_runs=3)
def test_rebuild_release(rebuilt_bin_path: Path):
    print(f"Rebuilding and installing last release to: {rebuilt_bin_path}")
    release = get_release(repository)
    assets = release["assets"]
    asset = next(iter([a for a in assets if a["name"] == get_asset_name(a)]), None)
    if not asset:
        warn(f"Couldn't find asset for OS {get_ostag()}, available assets:\n{assets}")

    with TemporaryDirectory() as td:
        # download the release
        download_path = Path(td)
        download_and_unzip(
            asset["browser_download_url"], path=download_path, verbose=True
        )

        # update IDEVELOPMODE
        source_files_path = download_path / asset["name"].replace(".zip", "") / "src"
        version_file_path = source_files_path / "Utilities" / "version.f90"
        with open(version_file_path) as f:
            lines = f.read().splitlines()
        assert len(lines) > 0, f"File is empty: {source_files_path}"
        with open(version_file_path, "w") as f:
            for line in lines:
                tag = "IDEVELOPMODE = 0"
                if tag in line:
                    line = line.replace(tag, "IDEVELOPMODE = 1")
                f.write(f"{line}\n")

        # rebuild with Meson
        meson_build(
            project_path=source_files_path.parent,
            build_path=download_path / "builddir",
            bin_path=rebuilt_bin_path,
        )


@flaky(max_runs=3)
def test_get_executables(downloaded_bin_path: Path):
    print(f"Installing MODFLOW-related executables to: {downloaded_bin_path}")
    downloaded_bin_path.mkdir(exist_ok=True, parents=True)
    flopy.utils.get_modflow(str(downloaded_bin_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Get executables needed for MODFLOW 6 testing")
    parser.add_argument("-p", "--path", help="path to top-level bin directory")
    args = parser.parse_args()
    bin_path = Path(args.path).resolve() if args.path else top_bin_path

    test_get_executables(bin_path / "downloaded")
    test_rebuild_release(bin_path / "rebuilt")

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory

import flopy
import pytest
from conftest import project_root_path
from flaky import flaky
from modflow_devtools.build import meson_build
from modflow_devtools.download import download_and_unzip, get_release
from modflow_devtools.misc import get_ostag

repository = "MODFLOW-ORG/modflow6"
top_bin_path = project_root_path / "bin"


def find_release_asset(release: dict) -> dict | None:
    """
    Find the platform-specific distribution archive for the current OS among
    a release's assets. Distribution archives are named "mf<version>_<ostag>.zip"
    (e.g. "mf6.8.0_linux.zip"). An exact-name match is required so that non-dist
    assets ("mf<version>_dfns.zip", "release.pdf") and other flavors
    ("mf<version>_win64ext.zip") are not selected. Intel macOS ("mac")
    distributions were dropped after 6.7.0, so fall back to the ARM build.
    """
    assets = release["assets"]
    version = release["tag_name"].lstrip("v")
    ostag = get_ostag()
    names = [f"mf{version}_{ostag}.zip"]
    if ostag == "mac":
        names.append(f"mf{version}_macarm.zip")
    for name in names:
        asset = next((a for a in assets if a["name"] == name), None)
        if asset is not None:
            return asset
    return None


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
    asset = find_release_asset(release)
    assert asset is not None, (
        f"Couldn't find a distribution asset for OS {get_ostag()}, available "
        f"assets:\n{[a['name'] for a in release['assets']]}"
    )

    with TemporaryDirectory() as td:
        # download the release
        download_path = Path(td)
        download_and_unzip(
            asset["browser_download_url"], path=download_path, verbose=True
        )

        # Force IDEVELOPMODE = 1 so DEV_* options used by "_dev" regression
        # models are accepted. Patch both version.f90 and version.f90.in: the
        # distribution's Meson build regenerates version.f90 from the template
        # via vcs_tag() at configure time, so patching version.f90 alone has
        # no effect.
        source_files_path = download_path / asset["name"].replace(".zip", "") / "src"
        version_dir = source_files_path / "Utilities"
        patched = False
        for version_file_path in (
            version_dir / "version.f90",
            version_dir / "version.f90.in",
        ):
            if not version_file_path.is_file():
                continue
            with open(version_file_path) as f:
                lines = f.read().splitlines()
            assert len(lines) > 0, f"File is empty: {version_file_path}"
            with open(version_file_path, "w") as f:
                for line in lines:
                    tag = "IDEVELOPMODE = 0"
                    if tag in line:
                        line = line.replace(tag, "IDEVELOPMODE = 1")
                        patched = True
                    f.write(f"{line}\n")
        assert patched, f"Failed to patch IDEVELOPMODE under {version_dir}"

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

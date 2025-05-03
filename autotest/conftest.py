import sys
from pathlib import Path
from warnings import warn

import pytest
from modflow_devtools.ostags import get_binary_suffixes

pytest_plugins = ["modflow_devtools.fixtures", "modflow_devtools.snapshots"]
project_root_path = Path(__file__).resolve().parent.parent


_exe_ext, _lib_ext = get_binary_suffixes(sys.platform)
_binaries_path = project_root_path / "bin"
_dl_bin_path = _binaries_path / "downloaded"
_rb_bin_path = _binaries_path / "rebuilt"
_binaries = {
    "development": [
        ("mf6", _binaries_path / f"mf6{_exe_ext}"),
        ("libmf6", _binaries_path / f"libmf6{_lib_ext}"),
        ("mf5to6", _binaries_path / f"mf5to6{_exe_ext}"),
        ("zbud6", _binaries_path / f"zbud6{_exe_ext}"),
    ],
    "downloaded": [
        ("mf2000", _dl_bin_path / f"mf2000{_exe_ext}"),
        ("mf2005", _dl_bin_path / f"mf2005dbl{_exe_ext}"),
        ("mfnwt", _dl_bin_path / f"mfnwtdbl{_exe_ext}"),
        ("mfusg", _dl_bin_path / f"mfusgdbl{_exe_ext}"),
        ("mflgr", _dl_bin_path / f"mflgrdbl{_exe_ext}"),
        ("mf2005s", _dl_bin_path / f"mf2005{_exe_ext}"),
        ("mt3dms", _dl_bin_path / f"mt3dms{_exe_ext}"),
        ("crt", _dl_bin_path / f"crt{_exe_ext}"),
        ("gridgen", _dl_bin_path / f"gridgen{_exe_ext}"),
        ("mp6", _dl_bin_path / f"mp6{_exe_ext}"),
        ("mp7", _dl_bin_path / f"mp7{_exe_ext}"),
        ("swtv4", _dl_bin_path / f"swtv4{_exe_ext}"),
        ("sutra", _dl_bin_path / f"sutra{_exe_ext}"),
        ("triangle", _dl_bin_path / f"triangle{_exe_ext}"),
        ("vs2dt", _dl_bin_path / f"vs2dt{_exe_ext}"),
        ("zonbudusg", _dl_bin_path / f"zonbudusg{_exe_ext}"),
    ],
    "rebuilt": [
        ("mf6_regression", _rb_bin_path / f"mf6{_exe_ext}"),
        ("libmf6_regression", _rb_bin_path / f"libmf6{_lib_ext}"),
        ("mf5to6_regression", _rb_bin_path / f"mf5to6{_exe_ext}"),
        ("zbud6_regression", _rb_bin_path / f"zbud6{_exe_ext}"),
    ],
}


@pytest.fixture(scope="session")
def bin_path() -> Path:
    """
    The directory containing binaries. Binaries under test live at the top
    level of this directory. Binaries downloaded from GitHub and/or built
    from the last official release are in subdirectories of this directory
    named "downloaded" and "rebuilt", respectively.
    """
    return _binaries_path


@pytest.fixture(scope="session")
def targets() -> dict[str, Path]:
    """
    Target executables for tests. These include local development builds as
    well as binaries 1) downloaded from GitHub and 2) rebuilt from the last
    official release.
    """

    d = dict()
    for k, v in _binaries["development"]:
        # require development binaries
        assert v.is_file(), f"Couldn't find binary '{k}' expected at: {v}"
        d[k] = v
    for k, v in _binaries["downloaded"]:
        # downloaded binaries are optional
        if v.is_file():
            d[k] = v
        else:
            warn(f"Couldn't find downloaded binary '{k}' expected at: {v}")
    for k, v in _binaries["rebuilt"]:
        # rebuilt binaries are optional
        if v.is_file():
            d[k] = v
        else:
            warn(f"Couldn't find rebuilt binary '{k}' expected at: {v}")
    return d


def try_get_target(targets: dict[str, Path], name: str) -> Path:
    """Try to retrieve the path to a binary. If the binary is a development
    target and can't be found, an error is raised. Otherwise (if the binary
    is downloaded or rebuilt) the test is skipped. This is to allow testing
    without downloaded or rebuilt binaries, e.g. if the network is down."""

    exe = targets.get(name)
    if exe:
        return exe
    elif name in _binaries["development"]:
        raise ValueError(f"Couldn't find binary '{name}'")
    else:
        pytest.skip(f"Couldn't find binary '{name}'")


@pytest.fixture
def original_regression(request) -> bool:
    """
    Whether to use the legacy regression test mechanism. This will run MF6
    models alongside non-MF6 model codes (e.g. MF2005) and compare results.

    The default regression test mechanism compares the MF6 under test with
    the latest MF6 release rebuilt in develop mode.
    """
    return request.config.getoption("--original-regression")


@pytest.fixture
def plot(request) -> bool:
    """
    Whether to make plots of model output. Useful for debugging.
    Disabled by default, enable with --plot on the command line.
    """
    return request.config.getoption("--plot")


@pytest.fixture(scope="session")
def markers(pytestconfig) -> str:
    """The markers used on the command line to filter tests."""
    return pytestconfig.getoption("-m")


@pytest.fixture(scope="session")
def models_path(request) -> list[Path]:
    """
    A directories containing model subdirectories. Use
    the --models-path command line option once or more to specify
    model directories. If at least one --models_path is provided,
    external tests (i.e. those using models from an external repo)
    will run against model input files found in the given location
    on the local filesystem rather than model input files from the
    official model registry. This is useful for testing changes to
    test model input files during MF6 development.
    """
    paths = request.config.getoption("--models-path") or []
    return [Path(p).expanduser().resolve().absolute() for p in paths]


def pytest_addoption(parser):
    parser.addoption(
        "--original-regression",
        action="store_true",
        default=False,
        help="use non-MF6 models for regression tests",
    )
    parser.addoption(
        "--parallel",
        action="store_true",
        default=False,
        help="include parallel test cases",
    )
    parser.addoption(
        "--netcdf",
        action="store_true",
        default=False,
        help="include netcdf test cases",
    )
    parser.addoption(
        "--plot",
        action="store_true",
        default=False,
        help="make plots of model output",
    )
    parser.addoption(
        "--models-path",
        action="append",
        type=str,
        help="directory containing model subdirectories. set this to run external "
        "tests (i.e. those using models from an external repo) against local model "
        "input files rather than input files from the official model registry.",
    )
    parser.addoption(
        "--namefile-pattern",
        action="store",
        type=str,
        default="mfsim.nam",
        help="namefile pattern to use when indexing models when --models-path is set."
        "does nothing otherwise. default is 'mfsim.nam'.",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--parallel"):
        skip_parallel = pytest.mark.skip(reason="need --parallel option to run")
        for item in items:
            if "parallel" in item.keywords:
                item.add_marker(skip_parallel)

    if not config.getoption("--netcdf"):
        skip_netcdf = pytest.mark.skip(reason="need --netcdf option to run")
        for item in items:
            if "netcdf" in item.keywords:
                item.add_marker(skip_netcdf)

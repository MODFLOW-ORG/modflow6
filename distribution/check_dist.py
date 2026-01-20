import platform
import subprocess
from os import environ
from pathlib import Path
from pprint import pprint

import pytest
from modflow_devtools.markers import no_parallel
from modflow_devtools.misc import run_cmd

# OS-specific extensions
SYSTEM = platform.system()
EXE_EXT = ".exe" if SYSTEM == "Windows" else ""
LIB_EXT = ".dll" if SYSTEM == "Windows" else ".so" if SYSTEM == "Linux" else ".dylib"
SCR_EXT = ".bat" if SYSTEM == "Windows" else ".sh"

# fortran compiler
FC = environ.get("FC", None)

# expected dirs
DIST_DIRS = {
    "develop": [
        "bin",
        "doc",
    ],
    "release": [
        "bin",
        "doc",
        "examples",
        "src",
        "srcbmi",
        "msvs",
        "make",
        "utils",
    ],
}


@pytest.fixture
def releasemode(request):
    return request.config.getoption("--releasemode")


@pytest.fixture
def dist_dir_path(request):
    def skip():
        pytest.skip(f"no distribution directory found at {path}")

    path = request.config.getoption("--path")
    if not path:
        skip()

    path = Path(path).expanduser().absolute()
    if not path.is_dir():
        skip()

    return path


@no_parallel
def test_directories(dist_dir_path, releasemode):
    for dir_path in DIST_DIRS["release" if releasemode else "develop"]:
        assert (dist_dir_path / dir_path).is_dir()


@no_parallel
def test_sources(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="sources not included in provisional distribution")

    # check top-level meson files
    assert (dist_dir_path / "meson.build").is_file()
    assert (dist_dir_path / "meson.options").is_file()

    # check src subdir
    assert (dist_dir_path / "src").is_dir()
    assert (dist_dir_path / "src" / "mf6.f90").is_file()
    version_file_path = dist_dir_path / "src" / "Utilities" / "version.f90"
    assert version_file_path.is_file()

    # check IDEVELOPMODE
    lines = open(version_file_path, "r").read().splitlines()
    pattern = ":: IDEVELOPMODE ="
    line = next(iter([l for l in lines if pattern in l]), None)
    assert line
    idevelopmode = 0 if releasemode else 1
    assert f"IDEVELOPMODE = {idevelopmode}" in line

    # check utils subdir
    assert (dist_dir_path / "utils").is_dir()
    assert (dist_dir_path / "utils" / "mf5to6").is_dir()
    assert (dist_dir_path / "utils" / "zonebudget").is_dir()
    assert (dist_dir_path / "utils" / "mf5to6" / "pymake").is_dir()
    assert (dist_dir_path / "utils" / "zonebudget" / "pymake").is_dir()
    assert not (dist_dir_path / "utils" / "idmloader").is_dir()


@no_parallel
@pytest.mark.skipif(not FC, reason="needs Fortran compiler")
def test_makefiles(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="makefiles not included in provisional distribution")

    assert (dist_dir_path / "make" / "makefile").is_file()
    assert (dist_dir_path / "make" / "makedefaults").is_file()
    assert (dist_dir_path / "utils" / "zonebudget" / "make" / "makefile").is_file()
    assert (dist_dir_path / "utils" / "zonebudget" / "make" / "makedefaults").is_file()
    assert (dist_dir_path / "utils" / "mf5to6" / "make" / "makefile").is_file()
    assert (dist_dir_path / "utils" / "mf5to6" / "make" / "makedefaults").is_file()

    # makefiles can't be used on Windows with ifort compiler
    if SYSTEM != "Windows" or FC != "ifort":
        print(subprocess.check_output("make", cwd=dist_dir_path / "make", shell=True))
        print(
            subprocess.check_output(
                "make",
                cwd=dist_dir_path / "utils" / "zonebudget" / "make",
                shell=True,
            )
        )
        print(
            subprocess.check_output(
                "make",
                cwd=dist_dir_path / "utils" / "mf5to6" / "make",
                shell=True,
            )
        )


@no_parallel
def test_msvs(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="MSVS files not included in provisional distribution")

    assert (dist_dir_path / "msvs" / "mf6.sln").is_file()
    assert (dist_dir_path / "msvs" / "mf6.vfproj").is_file()
    assert (dist_dir_path / "msvs" / "mf6bmi.sln").is_file()
    assert (dist_dir_path / "msvs" / "mf6bmi.vfproj").is_file()
    assert (dist_dir_path / "msvs" / "mf6core.vfproj").is_file()


@no_parallel
def test_docs(dist_dir_path, releasemode):
    # mf6io should always be included
    assert (dist_dir_path / "doc" / "mf6io.pdf").is_file()

    if releasemode:
        # check other custom-built documentation
        assert (dist_dir_path / "doc" / "release.pdf").is_file()
        assert (dist_dir_path / "doc" / "mf5to6.pdf").is_file()
        assert (dist_dir_path / "doc" / "zonebudget.pdf").is_file()
        assert (dist_dir_path / "doc" / "mf6suptechinfo.pdf").is_file()
        assert (dist_dir_path / "doc" / "mf6examples.pdf").is_file()

        # check publications downloaded from USGS site
        for pub in [
            "tm6a55",
            "tm6a56",
            "tm6a57",
            "tm6a61",
            "tm6a62",
        ]:
            assert (dist_dir_path / "doc" / f"{pub}.pdf").is_file()


@no_parallel
def test_examples(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="examples not included in provisional distribution")

    # check examples directory
    examples_path = dist_dir_path / "examples"
    assert examples_path.is_dir()

    # print examples found
    example_paths = [
        p for p in examples_path.glob("*") if p.is_dir() and p.stem.startswith("ex")
    ]
    assert any(example_paths)
    print(f"{len(example_paths)} example models found:")
    pprint(example_paths)

    # todo: check individual scripts? toggle via release workflow input?
    # model_paths = get_model_paths(examples_path)
    # script_paths = [mp / f"run{_scext}" for mp in model_paths]
    # for script_path in script_paths:
    #     print(f"Testing example script: {script_path}")
    #     assert script_path.is_file()
    #     out, err, ret = run_cmd(str(script_path), cwd=script_path.parent)
    #     assert not ret, out + err

    # check comprehensive examples script and give it a test run
    script_path = examples_path / f"runall{SCR_EXT}"
    print(f"Testing comprehensive examples script: {script_path}")
    assert script_path.is_file()
    out, err, ret = run_cmd(str(script_path), cwd=script_path.parent)
    assert not ret, out + err


@no_parallel
def test_binaries(dist_dir_path, releasemode):
    bin_path = dist_dir_path / "bin"
    assert (bin_path / f"mf6{EXE_EXT}").is_file()
    assert (bin_path / f"zbud6{EXE_EXT}").is_file()
    assert (bin_path / f"mf5to6{EXE_EXT}").is_file()
    assert (bin_path / f"libmf6{LIB_EXT}").is_file()

    # get version string
    output = " ".join(
        subprocess.check_output([str(bin_path / f"mf6{EXE_EXT}"), "-v"])
        .decode()
        .split()
    ).lower()
    assert output.startswith("mf6")

    # make sure version string reflects approval
    assert ("preliminary" in output) != releasemode

    # check version numbers
    version = output.lower().split(" ")[1]
    print("Version string:", version)
    v_split = version.split(".")
    assert len(v_split) >= 3

    # approved release should use semantic version number with
    # exactly 3 components and no alphabetic characters in it
    if releasemode:
        assert len(v_split) == 3
        assert all(s.isdigit() for s in v_split[:3])

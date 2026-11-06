import platform
import shutil
import subprocess
from os import environ
from pathlib import Path
from pprint import pprint

import pytest
from modflow_devtools.misc import get_model_paths

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
        "dfns",
        "doc",
    ],
    "release": [
        "bin",
        "dfns",
        "doc",
        "examples",
        "src",
        "srcbmi",
        "msvs",
        "make",
        "utils",
        "distribution",
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


def test_directories(dist_dir_path, releasemode):
    for dir_path in DIST_DIRS["release" if releasemode else "develop"]:
        assert (dist_dir_path / dir_path).is_dir()


def test_dfns(dist_dir_path):
    dfns_path = dist_dir_path / "dfns"
    assert dfns_path.is_dir()
    assert any(dfns_path.glob("*.dfn"))
    assert (dfns_path / "sim-nam.dfn").is_file()
    assert (dfns_path / "gwf-nam.dfn").is_file()


def test_sources(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="sources not included in provisional distribution")

    # check top-level meson files
    assert (dist_dir_path / "meson.build").is_file()
    assert (dist_dir_path / "meson.options").is_file()

    # check the version-tag helper invoked by the meson builds (issue #2973)
    assert (dist_dir_path / "distribution" / "vcs_tag_suffix.py").is_file()

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


@pytest.mark.xdist_group("build")
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


@pytest.mark.xdist_group("build")
@pytest.mark.skipif(not FC, reason="needs Fortran compiler")
def test_meson_build(dist_dir_path, releasemode):
    """Rebuild mf6 and mf5to6 from the archived sources with meson"""
    if not releasemode:
        pytest.skip(reason="sources not included in provisional distribution")

    # the helper the meson builds invoke for the version suffix
    assert (dist_dir_path / "distribution" / "vcs_tag_suffix.py").is_file()

    # base version string; in releasemode this is a bare "X.Y.Z"
    version_f90 = (dist_dir_path / "src" / "Utilities" / "version.f90").read_text()
    number_line = next(l for l in version_f90.splitlines() if ":: VERSIONNUMBER =" in l)
    base_version = number_line.split("'")[1]

    for label, setup_cwd, exe_rel in [
        ("mf6", dist_dir_path, Path("src") / f"mf6{EXE_EXT}"),
        (
            "mf5to6",
            dist_dir_path / "utils" / "mf5to6",
            Path("src") / f"mf5to6{EXE_EXT}",
        ),
    ]:
        build_dir = setup_cwd / "builddir"
        subprocess.check_call(["meson", "setup", "builddir"], cwd=setup_cwd)
        subprocess.check_call(["meson", "compile", "-C", "builddir"], cwd=setup_cwd)
        exe_path = build_dir / exe_rel
        assert exe_path.is_file(), f"{label} was not built at {exe_path}"

        if label == "mf6":
            # mf6 -v prints "<EXE NAME>: <VERSION>" (name includes the
            # extension on Windows) and exits 0
            tokens = subprocess.check_output([str(exe_path), "-v"]).decode().split()
            anchor = f"mf6{EXE_EXT}:"
        else:
            # mf5to6 has no version flag; its banner's second line is
            # "Version <FULLVERSION> [MM/DD/YYYY]". It then wants input, so
            # feed it EOF and ignore the (nonzero) exit status.
            proc = subprocess.run(
                [str(exe_path)],
                input="",
                capture_output=True,
                text=True,
                timeout=60,
                cwd=build_dir,
            )
            tokens = proc.stdout.split()
            anchor = "Version"
        print(f"{label} version string:", " ".join(tokens))

        # token after the anchor is the full version string
        reported = tokens[tokens.index(anchor) + 1]

        assert reported == base_version, (
            f"{label} reports {reported!r}, expected {base_version!r}"
        )


def test_msvs(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="MSVS files not included in provisional distribution")

    assert (dist_dir_path / "msvs" / "mf6.sln").is_file()
    assert (dist_dir_path / "msvs" / "mf6.vfproj").is_file()
    assert (dist_dir_path / "msvs" / "mf6bmi.sln").is_file()
    assert (dist_dir_path / "msvs" / "mf6bmi.vfproj").is_file()
    assert (dist_dir_path / "msvs" / "mf6core.vfproj").is_file()


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


def trace_script(script_path, mf6_path) -> list[Path]:
    """Example workspaces used in the run script, in order"""
    cwd = script_path.parent
    run_paths = []
    for line in script_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#!", "@echo", "echo", "pause")):
            continue
        if line.startswith("cd "):
            cwd = (cwd / line[3:].strip()).resolve()
            assert cwd.is_dir(), f"{script_path}: no such dir: {cwd}"
            continue
        exe_path = (cwd / line).resolve()
        assert exe_path == mf6_path.resolve(), f"{script_path}: not mf6: {line}"
        run_paths.append(cwd)
    return run_paths


def test_examples(dist_dir_path, releasemode):
    if not releasemode:
        pytest.skip(reason="examples not included in provisional distribution")

    # check examples directory
    examples_path = (dist_dir_path / "examples").resolve()
    assert examples_path.is_dir()

    # print examples found
    example_paths = [
        p for p in examples_path.glob("*") if p.is_dir() and p.stem.startswith("ex")
    ]
    assert any(example_paths)
    print(f"{len(example_paths)} example models found:")
    pprint(example_paths)

    # check the combined script runs all the examples
    mf6_path = dist_dir_path / "bin" / f"mf6{EXE_EXT}"
    runall_path = examples_path / f"runall{SCR_EXT}"
    assert runall_path.is_file()
    model_paths = trace_script(runall_path, mf6_path)
    assert len(model_paths) == len(set(model_paths))
    assert set(model_paths) == {p.resolve() for p in get_model_paths(examples_path)}

    # check each example's own script runs it
    for model_path in model_paths:
        script_path = model_path / f"run{SCR_EXT}"
        assert script_path.is_file()
        assert trace_script(script_path, mf6_path) == [model_path]


def get_examples(dist_path) -> dict[str, list[Path]]:
    """
    Map example names to their workspaces, in the order the distribution's
    combined run script runs them. Empty if the distribution has no examples
    or the script can't be traced; test_examples reports those problems.
    """
    examples_path = (dist_path / "examples").resolve()
    runall_path = examples_path / f"runall{SCR_EXT}"
    if not runall_path.is_file():
        return {}
    try:
        model_paths = trace_script(runall_path, dist_path / "bin" / f"mf6{EXE_EXT}")
    except AssertionError:
        return {}
    examples: dict[str, list[Path]] = {}
    for model_path in model_paths:
        name = model_path.relative_to(examples_path).parts[0]
        examples.setdefault(name, []).append(model_path)
    return examples


def pytest_generate_tests(metafunc):
    if "example_name" not in metafunc.fixturenames:
        return
    path = metafunc.config.getoption("--path")
    path = Path(path).expanduser().absolute() if path else None
    if metafunc.config.getoption("--releasemode") and path and path.is_dir():
        examples = get_examples(path)
    else:
        examples = {}
    if not examples:
        metafunc.parametrize(
            "example_name",
            [pytest.param(None, marks=pytest.mark.skip(reason="no examples found"))],
            ids=["none"],
        )
        return
    metafunc.parametrize("example_name", list(examples), ids=list(examples))


def test_example(dist_dir_path, releasemode, example_name, tmp_path):
    """
    Run an example with the distribution's mf6. The example is copied to a
    temp dir, keeping its workspaces' sibling layout so later workspaces can
    read earlier ones' output by relative path, and run in the given order.
    """
    if not releasemode:
        pytest.skip(reason="examples not included in provisional distribution")

    examples_path = (dist_dir_path / "examples").resolve()
    mf6_path = (dist_dir_path / "bin" / f"mf6{EXE_EXT}").absolute()
    shutil.copytree(examples_path / example_name, tmp_path / example_name)
    for model_path in get_examples(dist_dir_path)[example_name]:
        workspace = tmp_path / model_path.relative_to(examples_path)
        proc = subprocess.run(
            [mf6_path], cwd=workspace, capture_output=True, text=True, check=False
        )
        if proc.returncode or "Normal termination" not in proc.stdout:
            print(proc.stdout[-5000:], proc.stderr[-5000:])
            name = workspace.relative_to(tmp_path)
            pytest.fail(f"{name} failed (exit code {proc.returncode})")


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

    # make sure version string reflects build mode,
    # in develop mode it should have suffix "+sha"
    assert ("+" in output) != releasemode

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

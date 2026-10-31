"""
Run (and optionally benchmark) the MODFLOW 6 example models.

Unlike test_external_models.py, which tests each simulation workspace
in isolation, this module tests each *example*: all of an example's
workspaces are copied into the test's temporary directory, preserving
their sibling layout, and run in order. This matters for examples
where a later workspace reads an earlier one's output by relative
path (e.g. a GWT model's FMI package reading ../mf6gwf/flow.hds).

Each example runs with the development version (mf6) and with the
last release rebuilt in develop mode (mf6_regression), each in its
own copy of the example. Every workspace must terminate normally.

The tests only run when selected with -m examples, so they don't run
with the rest of the autotests, or with -m external and --models-path
pointing to other test models. Examples come from the model registry
by default, or from local directories with --models-path (e.g. the
examples/ folder in a distribution).

With --benchmark-enable, each example run is timed by pytest-benchmark,
grouped by example, with per-workspace times in the benchmark's extra
info. Use --benchmark-json to save the results, e.g. for the run-time
comparison table in the release notes.
"""

import subprocess
import time
from pathlib import Path

import pytest
from modflow_devtools.models import DEFAULT_REGISTRY, LocalRegistry
from modflow_devtools.ostags import get_binary_suffixes

# prefix of example models in the official model registry
PREFIX = "mf6/example/"
# examples to exclude
EXCLUDE = {
    # backward tracking examples read head and budget files reversed
    # by flopy after the flow model runs, which can't be done here.
    # they are also excluded from the distribution.
    "ex-prt-mp7-p02": "needs reversed head/budget files from flopy",
    "ex-prt-mp7-p04": "needs reversed head/budget files from flopy",
}
# mf6 binary for each version, relative to the bin path. this doesn't
# use the targets fixture, which requires binaries these tests don't.
EXE_EXT = get_binary_suffixes()[0]
VERSIONS = {
    "current": Path(f"mf6{EXE_EXT}"),
    "previous": Path("rebuilt", f"mf6{EXE_EXT}"),
}


def get_examples(registry) -> dict[str, list[str]]:
    """
    Map example names to the names of their workspaces in the registry,
    in the order they must run. Workspaces of a multi-model example are
    sorted with flow models first, so that transport, energy, particle
    tracking etc. models can consume flow model results. This is the
    same rule devtools uses to order examples in the official registry.
    """

    local = isinstance(registry, LocalRegistry)
    examples: dict[str, list[str]] = {}
    for model_name in registry.models:
        if not local and not model_name.startswith(PREFIX):
            continue
        example_name = model_name.removeprefix(PREFIX).split("/")[0]
        examples.setdefault(example_name, []).append(model_name)
    return {
        name: sorted(
            models,
            key=lambda m: 0 if "gwf" in m.removeprefix(PREFIX).partition("/")[2] else 1,
        )
        for name, models in sorted(examples.items())
    }


def pytest_generate_tests(metafunc):
    # see test_external_models.py and the --models-path option in conftest.py
    if "example_name" in metafunc.fixturenames:
        if "examples" not in (metafunc.config.getoption("-m") or ""):
            metafunc.parametrize("registry", [None], ids=["none"])
            metafunc.parametrize(
                "example_name",
                [pytest.param(None, marks=pytest.mark.skip(reason="need -m examples"))],
                ids=["none"],
            )
            return
        models_paths = metafunc.config.getoption("--models-path")
        models_paths = [
            Path(p).expanduser().resolve().absolute() for p in models_paths or []
        ]
        registry = LocalRegistry() if any(models_paths) else DEFAULT_REGISTRY
        registry_type = type(registry).__name__.lower().replace("registry", "")
        for path in models_paths:
            registry.index(path)
        examples = get_examples(registry)
        metafunc.parametrize("registry", [registry], ids=[registry_type])
        metafunc.parametrize("example_name", list(examples), ids=list(examples))


class RunError(Exception):
    pass


def run(exe: Path, workspace: Path) -> float:
    """Run the given mf6 binary in the workspace, return the elapsed time."""
    t0 = time.perf_counter()
    proc = subprocess.run(
        [exe], cwd=workspace, capture_output=True, text=True, check=False
    )
    elapsed = time.perf_counter() - t0
    if proc.returncode or "Normal termination" not in proc.stdout:
        print(proc.stdout[-5000:], proc.stderr[-5000:])
        raise RunError(
            f"{workspace.name} failed with {exe} (exit code {proc.returncode})"
        )
    return elapsed


@pytest.mark.slow
@pytest.mark.examples
@pytest.mark.parametrize("version", list(VERSIONS))
def test_example(registry, example_name, version, bin_path, function_tmpdir, benchmark):
    if example_name in EXCLUDE:
        pytest.skip(f"Skipping: {example_name} ({EXCLUDE[example_name]})")

    exe = (bin_path / VERSIONS[version]).absolute()
    if not exe.is_file():
        if version == "current":
            pytest.fail(f"Couldn't find binary: {exe}")
        pytest.skip(f"Couldn't find binary: {exe}")
    workspaces = []
    for model_name in get_examples(registry)[example_name]:
        workspace = function_tmpdir / model_name.removeprefix(PREFIX)
        registry.copy_to(workspace, model_name)
        workspaces.append(workspace)

    def run_example():
        return {str(ws.relative_to(function_tmpdir)): run(exe, ws) for ws in workspaces}

    benchmark.group = example_name
    benchmark.extra_info["version"] = version
    try:
        benchmark.extra_info["workspaces"] = benchmark.pedantic(
            run_example, rounds=1, iterations=1
        )
    except RunError as e:
        # the previous version is only a baseline for comparison, and
        # may not support everything the examples use. its failures
        # are reported as expected failures, not regressions.
        if version == "previous":
            pytest.xfail(str(e))
        pytest.fail(str(e))

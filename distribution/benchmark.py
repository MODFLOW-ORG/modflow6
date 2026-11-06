import argparse
import json
import subprocess
import sys
import textwrap
from os import PathLike
from pathlib import Path

from modflow_devtools.build import meson_build

from utils import get_project_root_path

PROJ_ROOT_PATH = get_project_root_path()
AUTOTEST_PATH = PROJ_ROOT_PATH / "autotest"
EXAMPLES_REPO_PATH = PROJ_ROOT_PATH.parent / "modflow6-examples"
BUILD_PATH = PROJ_ROOT_PATH / "builddir"
BIN_PATH = PROJ_ROOT_PATH / "bin"
BENCHMARKS_FILE_NAME = "run-time-comparison.md"
EXE_EXT = ".exe" if sys.platform.lower() == "win32" else ""


def get_mf6_info(exe: PathLike, flag: str) -> str:
    """Get the version (-v) or compiler (-c) string from an mf6 binary."""
    out = subprocess.run(
        [str(exe), flag], capture_output=True, text=True, check=True
    ).stdout.strip()
    return out.partition(":")[2].strip()


def elapsed_real_to_string(elt):
    if elt > 60.0:
        time_min = int(elt / 60.0)
        time_sec = elt % 60.0
        elt_str = f"{time_min} Minutes, "
    else:
        time_sec = elt
        elt_str = ""
    return elt_str + f"{time_sec:.3f} Seconds"


def get_binaries(build_path: PathLike, bin_path: PathLike):
    """
    Build the current version into the bin directory, and the last release
    into its "rebuilt" subdirectory, if they don't exist already.
    """

    bin_path = Path(bin_path).expanduser().absolute()
    if not (bin_path / f"mf6{EXE_EXT}").is_file():
        print("Building MODFLOW 6 development version")
        meson_build(
            project_path=PROJ_ROOT_PATH,
            build_path=Path(build_path).expanduser().absolute(),
            bin_path=bin_path,
        )
    if not (bin_path / "rebuilt" / f"mf6{EXE_EXT}").is_file():
        print("Rebuilding latest MODFLOW 6 release")
        subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "-v",
                "get_exes.py::test_rebuild_release",
                "--bin-path",
                str(bin_path),
            ],
            cwd=AUTOTEST_PATH,
            check=True,
        )


def run_benchmarks(
    bin_path: PathLike, examples_path: PathLike, json_path: PathLike
) -> bool:
    """
    Run the examples with the current and previous versions, benchmarking
    them with pytest-benchmark, and save results to the given JSON file.
    Return whether all examples ran successfully with the current version.
    """

    return not subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-v",
            "test_examples.py",
            "-m",
            "examples",
            "--models-path",
            str(Path(examples_path).expanduser().absolute()),
            "--bin-path",
            str(Path(bin_path).expanduser().absolute()),
            "--benchmark-enable",
            f"--benchmark-json={Path(json_path).expanduser().absolute()}",
        ],
        cwd=AUTOTEST_PATH,
        check=False,
    ).returncode


def write_results(bin_path: PathLike, json_path: PathLike, out_path: PathLike):
    """Write a markdown table of run times from the JSON benchmark results."""

    bin_path = Path(bin_path).expanduser().absolute()
    current_exe = bin_path / f"mf6{EXE_EXT}"
    previous_exe = bin_path / "rebuilt" / f"mf6{EXE_EXT}"
    current_v = get_mf6_info(current_exe, "-v").split()[0]
    previous_v = get_mf6_info(previous_exe, "-v").split()[0]

    # map workspace name to run time for each version
    times: dict[str, dict[str, float]] = {}
    benchmarks = json.loads(Path(json_path).read_text())["benchmarks"]
    for benchmark in benchmarks:
        version = benchmark["extra_info"]["version"]
        for workspace, t in benchmark["extra_info"]["workspaces"].items():
            # use the workspace name if it names the example scenario
            # (e.g. ex-gwf-hwell-a), else example/workspace (ex/mf6gwf)
            name = workspace.rpartition("/")[2]
            name = name if name.startswith("ex-") else workspace
            times.setdefault(name, {})[version] = t

    lines = [
        "### Benchmarks\n",
        "Comparison of run times of the current version of "
        f"MODFLOW 6 ({current_v}) to the previous version ({previous_v}). "
        "The current example models available from the "
        "[MODFLOW 6 Examples GitHub Repository]"
        "(https://github.com/MODFLOW-ORG/modflow6-examples) are "
        "used to compare run times. Simulations that fail are "
        "indicated by '--'. The percent difference, where calculated, "
        "is relative to the simulation run time for the previous "
        "version. Percent differences for example problems with "
        "short run times (less than 30 seconds) may not be significant.\n",
        f"{get_mf6_info(current_exe, '-c')}.\n\n",
        "| Example Problem "
        f"| Current Version {current_v} "
        f"| Previous Version {previous_v} "
        "| Percent difference |",
        "| :---------- | :----------: | :----------: | :----------: |",
    ]
    current_total = previous_total = 0.0
    for name, t in sorted(times.items()):
        current, previous = t.get("current"), t.get("previous")
        current_total += current or 0.0
        previous_total += previous or 0.0
        pd = (
            f"{(current - previous) / previous:.2%}"
            if current is not None and previous
            else "--"
        )
        lines.append(
            f"| {name} "
            f"| {'--' if current is None else elapsed_real_to_string(current)} "
            f"| {'--' if previous is None else elapsed_real_to_string(previous)} "
            f"| {pd} |"
        )
    pd = (current_total - previous_total) / previous_total
    lines.append(
        "| Total simulation time "
        f"| {elapsed_real_to_string(current_total)} "
        f"| {elapsed_real_to_string(previous_total)} "
        f"| {pd:.2%} |"
    )

    out_path = Path(out_path).expanduser().absolute()
    out_path.write_text("\n".join(lines) + "\n")
    print(f"Wrote benchmark results to {out_path}")


def benchmark(
    build_path: PathLike,
    bin_path: PathLike,
    examples_path: PathLike,
    out_path: PathLike,
):
    """
    Benchmark the current version of MODFLOW 6 against the last release
    with the example models, and write a markdown table of run times to
    the output directory. Raises an error if any example fails with the
    current version, after writing the table.
    """

    out_path = Path(out_path).expanduser().absolute()
    out_path.mkdir(parents=True, exist_ok=True)
    json_path = out_path / Path(BENCHMARKS_FILE_NAME).with_suffix(".json")
    get_binaries(build_path, bin_path)
    success = run_benchmarks(bin_path, examples_path, json_path)
    if not json_path.is_file():
        raise RuntimeError(f"No benchmark results, see pytest output: {json_path}")
    write_results(bin_path, json_path, out_path / BENCHMARKS_FILE_NAME)
    if not success:
        raise RuntimeError("Example(s) failed with the current version")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
    Benchmarks the current version of MODFLOW 6 against the latest official release,
    with the example models stored in the MODFLOW-ORG/modflow6-examples repository.
    The examples are run (and timed) by autotest/test_examples.py.
            """
        ),
    )
    parser.add_argument(
        "--build-path",
        required=False,
        default=str(BUILD_PATH),
        help="Path to the build workspace",
    )
    parser.add_argument(
        "--bin-path",
        required=False,
        default=str(BIN_PATH),
        help="Path to the directory containing current version executables, "
        "with the latest release in the 'rebuilt' subdirectory. Executables "
        "are built if not found.",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=False,
        default=str(PROJ_ROOT_PATH / "distribution"),
        help="Location to write the results",
    )
    parser.add_argument(
        "-e",
        "--examples-repo-path",
        required=False,
        default=str(EXAMPLES_REPO_PATH),
        help="Path to the directory with modflow6 examples",
    )
    args = parser.parse_args()
    examples_repo_path = Path(args.examples_repo_path)
    assert examples_repo_path.is_dir(), f"Examples repo not found: {examples_repo_path}"
    benchmark(
        build_path=args.build_path,
        bin_path=args.bin_path,
        examples_path=examples_repo_path / "examples",
        out_path=args.output_path,
    )

import argparse
import os
import shutil
import subprocess
import sys
import textwrap
from multiprocessing import Pool
from os import PathLike
from pathlib import Path

import flopy
import pytest
from modflow_devtools.build import meson_build
from modflow_devtools.download import download_and_unzip, get_latest_version
from modflow_devtools.misc import get_model_paths

from utils import get_project_root_path

PROJ_ROOT_PATH = get_project_root_path()
EXAMPLES_REPO_PATH = PROJ_ROOT_PATH.parent / "modflow6-examples"
BUILD_PATH = PROJ_ROOT_PATH / "builddir"
BIN_PATH = PROJ_ROOT_PATH / "bin"
GITHUB_REPO = "MODFLOW-ORG/modflow6"
BENCHMARKS_FILE_NAME = "run-time-comparison.md"
IS_WINDOWS = sys.platform.lower() == "win32"
EXE_EXT = ".exe" if IS_WINDOWS else ""
OSTAG = (
    "win64"
    if IS_WINDOWS
    else "linux"
    if sys.platform.lower().startswith("linux")
    else "mac"
)


def fetch_latest(outdir: PathLike) -> tuple[str, Path]:
    outdir = Path(outdir).expanduser().absolute()
    version = get_latest_version(GITHUB_REPO)
    distname = f"mf{version}_{OSTAG}"
    url = (
        f"https://github.com/{GITHUB_REPO}"
        + f"/releases/download/{version}/{distname}.zip"
    )
    download_and_unzip(url, path=outdir, verbose=True)
    return version, outdir / distname


def get_mf6_cmdargs(app, argv, text="mf6:", verbose=False):
    return_text = None
    proc = subprocess.Popen(
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=os.path.dirname(app)
    )
    result, _ = proc.communicate()
    if result is not None:
        c = result.decode("utf-8")
        c = c.rstrip("\r\n")
        if verbose:
            print(f"{c}")
        if text in c:
            idx0 = c.index(text) + len(text) + 1
            return_text = c[idx0:].strip()
    return return_text


def get_mf6_version(app, verbose=False):
    text = "mf6:"
    if str(app).endswith(".exe"):
        text = "mf6.exe:"
    version = get_mf6_cmdargs(app, [app, "-v"], text=text, verbose=verbose)
    if version is not None:
        version = version.split()[0]
        if verbose:
            print(f"version: {version}")
    return version


def get_mf6_compiler(app, verbose=False):
    text = "mf6:"
    if str(app).endswith(".exe"):
        text = "mf6.exe:"
    compiler = get_mf6_cmdargs(app, [app, "-c"], text=text, verbose=verbose)
    if verbose and compiler is not None:
        print(f"compiler: {compiler}")
    return compiler


def revert_files(app, example):
    """
    Modify input files to use deprecated/removed options
    as appropriate for the last release version
    """
    replacements = {".ims": {(6, 1, 1): ("dvclose", "hclose")}}
    extensions = list(replacements.keys())
    version = get_mf6_version(app)
    if version is not None:
        version = tuple([int(v) for v in version.split(".")])
    for file in os.listdir(example):
        _, extension = os.path.splitext(file)
        if extension in extensions:
            key = extension.lower()
            for v, replace in replacements[key].items():
                if version < v:
                    fpth = os.path.join(example, file)
                    with open(fpth, "r") as f:
                        lines = f.readlines()
                    with open(fpth, "w") as f:
                        for line in lines:
                            if replace[0] in line.lower():
                                line = line.lower().replace(replace[0], replace[1])
                            f.write(line)


def get_elapsed_time(buff, tag="Elapsed run time:"):
    elt_str = ""
    for line in buff:
        if tag in line:
            i0 = line.index(":")
            elt_str = line[i0 + 1 :].strip()
    return elt_str


def time_factor(time_unit):
    if "hours" in time_unit:
        factor = 60.0 * 60.0
    elif "minutes" in time_unit:
        factor = 60.0
    else:
        factor = 1.0
    return factor


def elapsed_string_to_real(elt_str):
    time_sec = 0.0
    t = elt_str.split()
    for idx in range(0, len(t), 2):
        t0 = float(t[idx])
        time_sec += t0 * time_factor(t[idx + 1].lower())
    return time_sec


def elapsed_real_to_string(elt):
    if elt > 60.0:
        time_min = int(elt / 60.0)
        time_sec = elt % 60.0
        elt_str = f"{time_min} Minutes, "
    else:
        time_sec = elt
        elt_str = ""
    return elt_str + f"{time_sec:.3f} Seconds"


def run_function(id, app, example):
    return (id, flopy.run_model(app, None, model_ws=example, silent=True, report=True))


def run_model(
    workspace: PathLike,
    dev_exe: PathLike,
    old_exe: PathLike,
):
    workspace = Path(workspace).expanduser().absolute()
    dev_exe = Path(dev_exe).expanduser().absolute()
    old_exe = Path(old_exe).expanduser().absolute()

    current_time = 0.0
    previous_time = 0.0
    generic_names = ["gwf", "gwt", "gwe", "prt"]
    generic_names = generic_names + [f"mf6{n}" for n in generic_names]
    name = (
        f"{workspace.parent.name}/{workspace.name}"
        if workspace.name in generic_names
        else workspace.name
    )
    print(f"Running scenario: {name}")
    line = f"| {name} |"
    prev_dir = os.path.join(workspace, "previous")
    if os.path.isdir(prev_dir):
        shutil.rmtree(prev_dir)
    print(f"Copying {workspace} ==> {prev_dir}")
    shutil.copytree(workspace, prev_dir)
    revert_files(old_exe, prev_dir)
    args = (
        (0, dev_exe, workspace),
        (1, old_exe, prev_dir),
    )
    pool = Pool(processes=2)
    results = [pool.apply_async(run_function, args=arg) for arg in args]
    pool.close()

    # set variables for processing
    id, (s, b) = results[0].get()
    if id == 0:
        success, buff = s, b
    elif id == 1:
        success0, buff0 = s, b
    id, (s, b) = results[1].get()
    if id == 0:
        success, buff = s, b
    elif id == 1:
        success0, buff0 = s, b

    if success:
        elt = get_elapsed_time(buff)
        line += f" {elt} |"
    else:
        print(f"Failure for current app with example: {name}")
        for b in buff:
            print(b)
        line += " -- |"

    if success0:
        elt0 = get_elapsed_time(buff0)
        line += f" {elt0} |"
    else:
        print(f"Failure for previous app with example: {name}")
        line += " -- |"

    if success and success0:
        t = elapsed_string_to_real(elt)
        t0 = elapsed_string_to_real(elt0)
        current_time += t
        previous_time += t0
        pd = (t - t0) / t0
        line += f" {pd:.2%} |"
    else:
        line += " -- |"

    # clean up previous directory
    if os.path.isdir(prev_dir):
        shutil.rmtree(prev_dir)

    return success, current_time, previous_time, line


def write_results(
    current_exe: PathLike,
    previous_exe: PathLike,
    output_path: PathLike,
    current_total,
    previous_total,
    lines: list[str],
):
    current_exe = Path(current_exe)
    previous_exe = Path(previous_exe)
    output_path = Path(output_path).expanduser().absolute()

    current_v = get_mf6_version(current_exe)
    previous_v = get_mf6_version(previous_exe)

    # open markdown table
    with open(output_path / BENCHMARKS_FILE_NAME, "w") as f:
        # get version numbers and write header

        line = "### Benchmarks\n\n"
        line += (
            "Comparison of run times of the current version of "
            + f"MODFLOW 6 ({current_v}) "
            + f"to the previous version ({previous_v}). "
            + "The current example models available from the "
            + "[MODFLOW 6 Examples GitHub Repository]"
            + "(https://github.com/MODFLOW-ORG/modflow6-examples) are "
            + "used to compare run times. Simulations that fail are "
            + "indicated by '--'. The percent difference, where calculated, "
            + "is relative to the simulation run time for the previous "
            + "version. Percent differences for example problems with "
            + "short run times (less than 30 seconds) may not be significant.\n\n"
            + f"{get_mf6_compiler(current_exe, verbose=True)}.\n\n\n"
        )
        line += "| Example Problem "
        line += f"| Current Version {current_v} "
        line += f"| Previous Version {previous_v} "
        line += "| Percent difference |\n"
        line += "| :---------- | :----------: | :----------: | :----------: |\n"
        f.write(line)

        # write benchmark data
        for line in lines:
            f.write(f"{line}\n")
            f.flush()

        # add final (total) line
        pd = (current_total - previous_total) / previous_total
        line = "| Total simulation time |"
        line += f" {elapsed_real_to_string(current_total)} |"
        line += f" {elapsed_real_to_string(previous_total)} |"
        line += f" {pd:.2%} |"
        f.write(f"{line}\n")


def run_benchmarks(
    build_path: PathLike,
    current_bin_path: PathLike,
    previous_bin_path: PathLike,
    examples_path: PathLike,
    out_path: PathLike,
    excluded: list[str] = [],
):
    """Benchmark current development version against previous release
    with example models.
    """

    build_path = Path(build_path).expanduser().absolute()
    current_bin_path = Path(current_bin_path).expanduser().absolute()
    previous_bin_path = Path(previous_bin_path).expanduser().absolute()
    examples_path = Path(examples_path).expanduser().absolute()
    out_path = Path(out_path).expanduser().absolute()

    example_dirs = get_model_paths(examples_path, excluded=excluded)
    if not any(example_dirs):
        raise FileNotFoundError(f"No example model paths found in {examples_path}")

    exe_name = f"mf6{EXE_EXT}"
    dev_exe = current_bin_path / exe_name
    old_exe = previous_bin_path / exe_name

    if not dev_exe.is_file():
        print("Building MODFLOW 6 development version")
        meson_build(
            project_path=PROJ_ROOT_PATH,
            build_path=build_path,
            bin_path=current_bin_path,
        )

    if not old_exe.is_file():
        version, download_path = fetch_latest(out_path)
        print(f"Rebuilding latest MODFLOW 6 release {version} in development mode")
        meson_build(
            project_path=download_path,
            build_path=build_path,
            bin_path=previous_bin_path,
        )
        shutil.rmtree(download_path)

    print("Benchmarking MODFLOW 6 versions:")
    print(f"    dev: {dev_exe}")
    print(f"    old: {old_exe}")

    # benchmark models
    current_total = 0.0
    previous_total = 0.0
    lines = []
    skip = ["ex-prt-mp7-p02", "ex-prt-mp7-p04"]
    for example_dir in example_dirs:
        if any(
            (pattern in example_dir.name or pattern in example_dir.parent.name)
            for pattern in skip
        ):
            print(f"Skipping {example_dir}")
            continue
        success, t, t0, line = run_model(
            workspace=example_dir, dev_exe=dev_exe, old_exe=old_exe
        )
        if not success:
            print(f"{example_dir} run failed")
        current_total += t
        previous_total += t0
        lines.append(line)

    # create markdown results file
    write_results(
        current_exe=dev_exe,
        previous_exe=old_exe,
        output_path=out_path,
        current_total=current_total,
        previous_total=previous_total,
        lines=lines,
    )


@pytest.mark.skip(reason="for manual testing")
def test_run_benchmarks(tmp_path):
    run_benchmarks(
        build_path=BUILD_PATH,
        current_bin_path=BIN_PATH,
        previous_bin_path=BIN_PATH / "rebuilt",
        examples_path=EXAMPLES_REPO_PATH / "examples",
        out_path=tmp_path,
        excluded=["previous"],
    )
    assert (tmp_path / BENCHMARKS_FILE_NAME).is_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
    Benchmarks the current version of MODFLOW 6 against the latest official release,
    with the example models stored in the MODFLOW-ORG/modflow6-examples repository.
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
        "--current-bin-path",
        required=False,
        default=str(BIN_PATH),
        help="Path to the directory to install current version executables",
    )
    parser.add_argument(
        "--previous-bin-path",
        required=False,
        default=str(BIN_PATH / "rebuilt"),
        help="Path to the directory to install previous version executables",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=False,
        default=str(PROJ_ROOT_PATH / "distribution" / ""),
        help="Location to create the zip archive",
    )
    parser.add_argument(
        "-e",
        "--examples-repo-path",
        required=False,
        default=str(PROJ_ROOT_PATH.parent / "modflow6-examples"),
        help="Path to the directory with modflow6 examples",
    )
    args = parser.parse_args()
    build_path = Path(args.build_path)
    current_bin_path = Path(args.current_bin_path)
    previous_bin_path = Path(args.previous_bin_path)
    output_path = Path(args.output_path) if args.output_path else Path(os.getcwd())
    examples_repo_path = (
        Path(args.examples_repo_path) if args.examples_repo_path else EXAMPLES_REPO_PATH
    )

    output_path.mkdir(parents=True, exist_ok=True)
    assert examples_repo_path.is_dir(), f"Examples repo not found: {examples_repo_path}"

    run_benchmarks(
        build_path=build_path,
        current_bin_path=current_bin_path,
        previous_bin_path=previous_bin_path,
        examples_path=examples_repo_path / "examples",
        out_path=output_path,
        excluded=["previous"],
    )

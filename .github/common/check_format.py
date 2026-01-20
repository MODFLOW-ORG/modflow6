import argparse
import os
import sys
import timeit
from itertools import repeat
from multiprocessing import Manager, Pool, cpu_count
from pathlib import Path
from queue import Empty
from subprocess import run

PROJ_ROOT = Path(__file__).parents[2]

# exclude these directories from checks
excludedirs = [
    PROJ_ROOT / ".pixi",
    PROJ_ROOT / "src" / "Utilities" / "Libraries" / "blas",
    PROJ_ROOT / "src" / "Utilities" / "Libraries" / "daglib",
    PROJ_ROOT / "src" / "Utilities" / "Libraries" / "rcm",
    PROJ_ROOT / "src" / "Utilities" / "Libraries" / "sparsekit",
    PROJ_ROOT / "src" / "Utilities" / "Libraries" / "sparskit2",
    PROJ_ROOT / "utils" / "mf5to6",
    PROJ_ROOT / "utils" / "zonebudget" / "msvs",
    PROJ_ROOT / "msvs",
]

# exclude these files from checks
excludefiles = [
    PROJ_ROOT / "src" / "Idm" / "gwf-csubidm.f90",
    PROJ_ROOT / "src" / "Idm" / "gwf-stoidm.f90",
    PROJ_ROOT / "src" / "Idm" / "gwf-vscidm.f90",
]

# commands
fprettify = "fprettify -c .fprettify.yaml"


def excluded(path) -> bool:
    path = Path(path)
    for f in excludefiles:
        if os.path.exists(f) and os.path.samefile(path, f):
            return True
    for d in excludedirs:
        if os.path.exists(d) and path.is_relative_to(d):
            return True
    return False


def check_format(path, lock, checks, failures, write_changes=False, verbose=False):
    path = Path(path)
    if excluded(path):
        if verbose:
            print(f"Skipping format check: {path}")
        return

    if verbose:
        print(f"Checking format: {path}")

    diff = "" if write_changes else "-d"
    cmd = f"{fprettify} {diff} {path}"
    result = run(cmd, capture_output=True, shell=True)
    if result.stdout or result.stderr:
        failures.put(path)

    with lock:
        checks.value += 1


def report(checks, failures, duration: float) -> bool:
    def pop(q):
        return q.get(block=False)

    n_failures = failures.qsize()
    success = n_failures == 0
    print(f"Checked format for {checks.value} Fortran files in {duration:.4f}s")

    if n_failures > 0:
        success = False
        stats = f"failures: {n_failures} \\"
        hr = "".join(repeat("_", len(stats) - 1))
        print(f"{hr}\n{stats}")
        while True:
            try:
                print(f"{fprettify} {pop(failures)}")
            except Empty:
                break

    print()
    return success


if __name__ == "__main__":
    start = timeit.default_timer()
    parser = argparse.ArgumentParser(
        "MODFLOW 6 Fortran source code format verification"
    )
    parser.add_argument(
        "-p",
        "--path",
        help="path to file or directory",
        default=PROJ_ROOT,
    )
    parser.add_argument(
        "-w",
        "--write-changes",
        help="write codespell changes in place if possible",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="verbose", default=False
    )
    args = parser.parse_args()
    path = Path(args.path).expanduser().absolute()
    assert path.exists(), f"Path not found: {path}"
    write = args.write_changes
    verbose = args.verbose

    # shared state
    manager = Manager()
    lock = manager.Lock()
    checks = manager.Value("checks", 0)
    failures = manager.Queue()

    if path.is_file():
        check_format(path, lock, checks, failures, verbose=verbose)
    else:
        with Pool(cpu_count()) as pool:
            files = [str(p) for p in path.rglob("*.[fF]9[05]")]
            if verbose:
                msg = f"Checking {len(files)} files in directory: {path}"
                print(msg)
                print("".join(repeat("-", len(msg))))
            pool.starmap(
                check_format,
                [(f, lock, checks, failures, write, verbose) for f in files],
            )

    stop = timeit.default_timer()
    sys.exit(0 if report(checks, failures, stop - start) else 1)

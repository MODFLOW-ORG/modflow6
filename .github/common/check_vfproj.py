import argparse
import xml.etree.ElementTree as ET
from itertools import chain
from pathlib import Path
from pprint import pformat

PROJ_ROOT = Path(__file__).parents[2]


def get_source_files(src_path, verbose=False):
    if verbose:
        print(f"    Checking source files in dir: {src_path}")
    extensions = ("*.[fF]", "*.[fF]9[05]", "*.inc")
    for ext in extensions:
        for path in src_path.glob(f"**/{ext}"):
            yield path.absolute()


def get_extra_files(extrafiles_path, src_path, extra_path, verbose=False):
    if verbose:
        print(f"    Checking extra files in dir: {extra_path}")
    paths = set(get_source_files(extra_path))
    with open(extrafiles_path) as ef:
        for path in ef:
            path = path.replace("\\", "/").strip()
            root = extra_path if "../../../" in path else src_path
            path = root / path.replace("../", "").replace("src/", "")
            if path in paths:
                yield path.absolute()


def get_msvs_files(vfproj_path, src_path, extra_path=None, verbose=False):
    if verbose:
        print(f"    Checking MSVS file: {vfproj_path}")
    tree = ET.parse(vfproj_path)
    root = tree.getroot()
    for f in root.iter("File"):
        path = f.attrib["RelativePath"].replace("\\", "/")
        yield (
            (
                extra_path
                if (extra_path is not None and "../../../" in path)
                else src_path
            )
            / path.replace("../", "").replace("src/", "").replace("srcbmi/", "")
        ).absolute()


def check_files(name, src, msvs, verbose=False):
    if verbose:
        print(f"Checking that {name} MSVS files match source files...")
    s, m = set(src), set(msvs)
    diff = s ^ m
    assert not any(diff), (
        f"{name} src files don't match MSVS project file\n"
        f"=> symmetric difference:\n{pformat(diff)}\n"
        f"=> src - msvs:\n{pformat(s - m)}\n"
        f"=> msvs - src:\n{pformat(m - s)}\n"
    )


def check_mf6(verbose):
    src_path = PROJ_ROOT / "src"
    src = get_source_files(src_path=src_path, verbose=verbose)
    msvs = chain(
        *[
            get_msvs_files(
                vfproj_path=PROJ_ROOT / "msvs" / f, src_path=src_path, verbose=verbose
            )
            for f in ["mf6core.vfproj", "mf6.vfproj"]
        ]
    )
    check_files("MF6", src, msvs, verbose)


def check_bmi(verbose):
    src_path = PROJ_ROOT / "srcbmi"
    src = get_source_files(src_path=src_path, verbose=verbose)
    msvs = get_msvs_files(
        vfproj_path=PROJ_ROOT / "msvs" / "mf6bmi.vfproj",
        src_path=src_path,
        verbose=verbose,
    )
    check_files("BMI", src, msvs, verbose)


def check_mf5to6(verbose):
    util_path = PROJ_ROOT / "utils" / "mf5to6"
    src_path = util_path / "src"
    src_files = chain(
        get_source_files(src_path=src_path, verbose=verbose),
        get_extra_files(
            extrafiles_path=util_path / "pymake" / "extrafiles.txt",
            src_path=src_path,
            extra_path=PROJ_ROOT / "src",
            verbose=verbose,
        ),
    )
    msvs_files = get_msvs_files(
        vfproj_path=util_path / "msvs" / "mf5to6.vfproj",
        src_path=src_path,
        extra_path=PROJ_ROOT / "src",
        verbose=verbose,
    )
    check_files("MODFLOW 5 to 6 converter", src_files, msvs_files, verbose)


def check_zonebudget(verbose):
    util_path = PROJ_ROOT / "utils" / "zonebudget"
    src_path = util_path / "src"
    src_files = chain(
        get_source_files(src_path=src_path, verbose=verbose),
        get_extra_files(
            extrafiles_path=util_path / "pymake" / "extrafiles.txt",
            src_path=src_path,
            extra_path=PROJ_ROOT / "src",
            verbose=verbose,
        ),
    )
    msvs_files = get_msvs_files(
        vfproj_path=util_path / "msvs" / "zonebudget.vfproj",
        src_path=src_path,
        extra_path=PROJ_ROOT / "src",
        verbose=verbose,
    )
    check_files("Zonebudget", src_files, msvs_files, verbose)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Check Microsoft Visual Studio project files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=True,
        required=False,
        help="Show verbose output",
    )
    args = parser.parse_args()
    verbose = args.verbose
    check_mf6(verbose)
    check_bmi(verbose)
    check_mf5to6(verbose)
    check_zonebudget(verbose)
    print("MSVS project (vfproj) files are up-to-date.")

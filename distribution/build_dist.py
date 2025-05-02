import argparse
import os
import platform
import textwrap
from os import PathLike, environ
from pathlib import Path
from pprint import pprint
from shutil import copy, copyfile, copytree, ignore_patterns, rmtree
from typing import Optional

import pytest
from build_docs import build_documentation
from build_makefiles import (
    build_mf5to6_makefile,
    build_mf6_makefile,
    build_zbud6_makefile,
)
from modflow_devtools.build import meson_build
from modflow_devtools.download import download_and_unzip, get_release
from modflow_devtools.markers import no_parallel, requires_exe
from modflow_devtools.misc import get_model_paths

from utils import get_project_root_path

# default paths
PROJ_ROOT_PATH = get_project_root_path()
BUILDDIR_PATH = PROJ_ROOT_PATH / "builddir"
DEFAULT_MODELS = ["gwf", "gwt", "gwe", "prt", "swf"]

# OS-specific extensions
SYSTEM = platform.system()
EXE_EXT = ".exe" if SYSTEM == "Windows" else ""
LIB_EXT = ".dll" if SYSTEM == "Windows" else ".so" if SYSTEM == "Linux" else ".dylib"
SCR_EXT = ".bat" if SYSTEM == "Windows" else ".sh"
MF6_EXE = f"mf6{EXE_EXT}"

# Fortran and C compilers
FC = environ.get("FC", "gfortran")
CC = environ.get("CC", "gcc")


def copy_sources(output_path: PathLike):
    output_path = Path(output_path).expanduser().absolute()

    # make sure output directory exists
    output_path.mkdir(exist_ok=True)

    # Copy Visual Studio sln and project files
    print("Copying msvs files to output directory")
    (output_path / "msvs").mkdir(exist_ok=True)
    source_msvs_path = PROJ_ROOT_PATH / "msvs"
    for d in [
        str(source_msvs_path / "mf6.sln"),
        str(source_msvs_path / "mf6.vfproj"),
        str(source_msvs_path / "mf6core.vfproj"),
        str(source_msvs_path / "mf6bmi.sln"),
        str(source_msvs_path / "mf6bmi.vfproj"),
    ]:
        copy(d, output_path / "msvs")

    ignored = [".DS_Store"]

    # copy top-level meson.build and meson.options
    copy(PROJ_ROOT_PATH / "meson.build", output_path)
    copy(PROJ_ROOT_PATH / "meson.options", output_path)

    # copy source folder
    src_path = PROJ_ROOT_PATH / "src"
    dst_path = output_path / "src"
    print(f"Copying {src_path} to {dst_path}")
    copytree(src_path, dst_path, ignore=ignore_patterns(*ignored))

    # copy srcbmi folder
    src_path = PROJ_ROOT_PATH / "srcbmi"
    dst_path = output_path / "srcbmi"
    print(f"Copying {src_path} to {dst_path}")
    copytree(src_path, dst_path, ignore=ignore_patterns(*ignored))

    # copy utils folder
    src_path = PROJ_ROOT_PATH / "utils"
    dst_path = output_path / "utils"
    print(f"Copying {src_path} to {dst_path}")
    ignored.extend(["idmloader"])
    copytree(src_path, dst_path, ignore=ignore_patterns(*ignored))


@no_parallel
def test_copy_sources(tmp_path):
    copy_sources(tmp_path)

    assert (tmp_path / "src").is_dir()
    assert (tmp_path / "srcbmi").is_dir()
    assert (tmp_path / "utils").is_dir()
    assert (tmp_path / "msvs").is_dir()

    assert (tmp_path / "src" / "meson.build").is_file()
    assert (tmp_path / "srcbmi" / "meson.build").is_file()
    assert (tmp_path / "utils" / "mf5to6" / "meson.build").is_file()
    assert (tmp_path / "msvs" / "mf6.sln").is_file()

    assert (tmp_path / "utils").is_dir()
    assert (tmp_path / "utils" / "mf5to6").is_dir()
    assert (tmp_path / "utils" / "zonebudget").is_dir()
    assert (tmp_path / "utils" / "mf5to6" / "pymake").is_dir()
    assert (tmp_path / "utils" / "zonebudget" / "pymake").is_dir()
    assert not (tmp_path / "utils" / "idmloader").is_dir()


def setup_examples(
    bin_path: PathLike,
    examples_path: PathLike,
    force: bool = False,
    models: Optional[list[str]] = None,
):
    examples_path = Path(examples_path).expanduser().absolute()

    # find and download example models distribution from latest examples release
    latest = get_release("MODFLOW-ORG/modflow6-examples", tag="latest", verbose=True)
    assets = latest["assets"]
    print(f"Found {len(assets)} assets from the latest examples release:")
    pprint([a["name"] for a in assets])
    asset = next(iter([a for a in assets if a["name"].endswith("examples.zip")]), None)
    download_and_unzip(asset["browser_download_url"], examples_path, verbose=True)

    # filter examples for models selected for release
    # and omit any excluded models
    excluded = ["ex-prt-mp7-p02", "ex-prt-mp7-p04"]
    for p in examples_path.glob("*"):
        if not any(m in p.stem for m in models):
            print(f"Omitting example due to model selection: {p.stem}")
            rmtree(p)
        if any(e in p.stem for e in excluded):
            print(f"Omitting deliberately excluded example: {p.stem}")
            rmtree(p)

    # list folders with mfsim.nam (recursively)
    # and add run.sh/bat script to each folder
    model_paths = get_model_paths(examples_path)
    for mp in model_paths:
        script_path = mp / f"run{SCR_EXT}"
        if not force and script_path.is_file():
            print(f"Script {script_path} already exists")
        else:
            print(f"Creating {script_path}")
            with open(script_path, "w") as f:
                if SYSTEM == "Windows":
                    f.write("@echo off" + "\n")
                else:
                    f.write("#!/bin/sh" + "\n")
                runbatloc = os.path.relpath(bin_path / MF6_EXE, start=mp)
                f.write(runbatloc + "\n")
                if SYSTEM == "Windows":
                    f.write("echo." + "\n")
                    f.write("echo Run complete.  Press any key to continue" + "\n")
                    f.write("pause>nul" + "\n")

            if SYSTEM != "Windows":
                script_path.chmod(script_path.stat().st_mode | 0o111)
                print(f"Execute permission set for {script_path}")

    # add runall.sh/bat, which runs all examples
    script_path = examples_path / f"runall{SCR_EXT}"
    if not force and script_path.is_file():
        print(f"Script {script_path} already exists")
    else:
        print(f"Creating {script_path}")
        with open(script_path, "w") as f:
            if SYSTEM != "Windows":
                f.write("#!/bin/sh" + "\n")
            for mp in model_paths:
                d = os.path.relpath(mp, start=examples_path)
                s = f"cd {d}"
                f.write(s + "\n")
                runbatloc = os.path.relpath(bin_path / MF6_EXE, start=mp)
                f.write(runbatloc + "\n")
                d = os.path.relpath(examples_path, start=mp)
                s = f"cd {d}"
                f.write(s + "\n")
                s = ""
                f.write(s + "\n")
            if SYSTEM == "Windows":
                f.write("pause" + "\n")
            else:
                script_path.chmod(script_path.stat().st_mode | 0o111)
                print(f"Execute permission set for {script_path}")


def build_programs_meson(build_path: PathLike, bin_path: PathLike, force: bool = False):
    build_path = Path(build_path).expanduser().absolute()
    bin_path = Path(bin_path).expanduser().absolute()

    exe_paths = [
        bin_path / f"mf6{EXE_EXT}",
        bin_path / f"zbud6{EXE_EXT}",
        bin_path / f"mf5to6{EXE_EXT}",
    ]
    lib_paths = [bin_path / f"libmf6{LIB_EXT}"]

    if (
        not force
        and all(p.is_file() for p in exe_paths)
        and all(p.is_file() for p in lib_paths)
    ):
        print("Binaries already exist:")
        pprint(exe_paths + lib_paths)
    else:
        print(f"Building binaries in {build_path}, installing to {bin_path}")
        meson_build(
            project_path=PROJ_ROOT_PATH,
            build_path=build_path,
            bin_path=bin_path,
        )
        meson_build(
            project_path=PROJ_ROOT_PATH / "utils/mf5to6",
            build_path=build_path,
            bin_path=bin_path,
        )

    for target in exe_paths + lib_paths:
        assert target.is_file(), f"Failed to build {target}"
        target.chmod(target.stat().st_mode | 0o111)
        print(f"Execute permission set for {target}")


@no_parallel
def test_build_programs_meson(tmp_path):
    build_programs_meson(tmp_path / "builddir", tmp_path / "bin")


def build_makefiles(output_path: PathLike):
    output_path = Path(output_path).expanduser().absolute()

    # create and copy mf6 makefile
    build_mf6_makefile()
    (output_path / "make").mkdir(parents=True, exist_ok=True)
    copyfile(
        PROJ_ROOT_PATH / "make" / "makefile",
        output_path / "make" / "makefile",
    )
    copyfile(
        PROJ_ROOT_PATH / "make" / "makedefaults",
        output_path / "make" / "makedefaults",
    )

    # create and copy zbud6 makefile
    build_zbud6_makefile()
    rel_path = Path("utils") / "zonebudget" / "make"
    (output_path / rel_path).mkdir(parents=True, exist_ok=True)
    copyfile(
        PROJ_ROOT_PATH / rel_path / "makefile",
        output_path / rel_path / "makefile",
    )
    copyfile(
        PROJ_ROOT_PATH / rel_path / "makedefaults",
        output_path / rel_path / "makedefaults",
    )

    # create and copy mf5to6 makefile
    build_mf5to6_makefile()
    rel_path = Path("utils") / "mf5to6" / "make"
    (output_path / rel_path).mkdir(parents=True, exist_ok=True)
    copyfile(
        PROJ_ROOT_PATH / rel_path / "makefile",
        output_path / rel_path / "makefile",
    )
    copyfile(
        PROJ_ROOT_PATH / rel_path / "makedefaults",
        output_path / rel_path / "makedefaults",
    )


@no_parallel
def test_build_makefiles(tmp_path):
    build_makefiles(tmp_path)

    assert (tmp_path / "make" / "makefile").is_file()
    assert (tmp_path / "make" / "makedefaults").is_file()
    assert (tmp_path / "utils" / "zonebudget" / "make" / "makefile").is_file()
    assert (tmp_path / "utils" / "zonebudget" / "make" / "makedefaults").is_file()
    assert (tmp_path / "utils" / "mf5to6" / "make" / "makefile").is_file()
    assert (tmp_path / "utils" / "mf5to6" / "make" / "makedefaults").is_file()

    os.system(f"cd {tmp_path} && make -f make/makefile")
    os.system(f"cd {tmp_path} && make -f utils/zonebudget/make/makefile")
    os.system(f"cd {tmp_path} && make -f utils/mf5to6/make/makefile")


def build_distribution(
    build_path: PathLike,
    output_path: PathLike,
    full: bool = False,
    force: bool = False,
    models: Optional[list[str]] = None,
):
    print(f"Building {'full' if full else 'minimal'} distribution")

    build_path = Path(build_path).expanduser().absolute()
    output_path = Path(output_path).expanduser().absolute()

    # binaries
    build_programs_meson(
        build_path=build_path,
        bin_path=output_path / "bin",
        force=force,
    )

    # code.json metadata
    copy(PROJ_ROOT_PATH / "code.json", output_path)

    # full releases include examples, source code, makefiles and docs
    if not full:
        return

    # download and setup example models
    setup_examples(
        bin_path=output_path / "bin",
        examples_path=output_path / "examples",
        force=force,
        models=models,
    )

    # copy source code files
    copy_sources(output_path=output_path)

    # build and copy makefiles
    build_makefiles(output_path=output_path)

    # build docs
    build_documentation(
        bin_path=output_path / "bin",
        full=full,
        out_path=output_path / "doc",
        force=force,
    )


@no_parallel
@requires_exe("pdflatex")
@pytest.mark.skip(reason="manual testing")
@pytest.mark.parametrize("full", [True, False])
def test_build_distribution(tmp_path, full):
    output_path = tmp_path / "dist"
    build_distribution(
        build_path=tmp_path / "builddir",
        output_path=output_path,
        full=full,
        force=True,
    )

    if full:
        # todo
        pass
    else:
        # check binaries and libs
        system = platform.system()
        ext = ".exe" if system == "Windows" else ""
        for exe in ["mf6", "mf5to6", "zbud6"]:
            assert (output_path / f"{exe}{ext}").is_file()
        assert (
            output_path
            / (
                "libmf6"
                + (
                    ".so"
                    if system == "Linux"
                    else (".dylib" if system == "Darwin" else ".dll")
                )
            )
        ).is_file()

        # check mf6io docs
        assert (output_path / "mf6io.pdf").is_file()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Create a MODFLOW 6 distribution. If output path is provided
            distribution files are written to the selected path, if not
            they are written to the distribution/ project subdirectory.
            By default a minimal distribution containing only binaries,
            mf6io documentation, release notes and metadata (code.json)
            is created. To create a full distribution including sources
            and examples, use the --full flag. Models to be included in
            the examples and documentation can be selected with --model
            (or -m), which may be used multiple times. Use --force (-f)
            to overwrite preexisting distribution artifacts; by default
            the script is lazy and will only create what it can't find.
            """
        ),
    )
    parser.add_argument(
        "--build-path",
        required=False,
        default=str(BUILDDIR_PATH),
        help="Path to the build workspace",
    )
    parser.add_argument(
        "-o",
        "--output-path",
        required=False,
        default=str(PROJ_ROOT_PATH / "distribution"),
        help="Path to create distribution artifacts",
    )
    parser.add_argument(
        "-m",
        "--model",
        required=False,
        action="append",
        help="Filter models to include",
    )
    parser.add_argument(
        "--full",
        required=False,
        default=False,
        action="store_true",
        help="Build a full rather than minimal distribution",
    )
    parser.add_argument(
        "-f",
        "--force",
        required=False,
        default=False,
        action="store_true",
        help="Recreate and overwrite existing artifacts",
    )
    args = parser.parse_args()
    build_path = Path(args.build_path)
    out_path = Path(args.output_path)
    out_path.mkdir(parents=True, exist_ok=True)
    models = args.model if args.model else DEFAULT_MODELS

    build_distribution(
        build_path=build_path,
        output_path=out_path,
        full=args.full,
        force=args.force,
        models=models,
    )

import os
import sys
from os import environ
from pathlib import Path

import pytest
from flaky import flaky
from modflow_devtools.markers import no_parallel, requires_exe
from modflow_devtools.misc import set_dir

import pymake
from utils import get_modified_time, get_project_root_path

PROJ_ROOT_PATH = get_project_root_path()
IS_WINDOWS = sys.platform.lower() == "win32"
EXE_EXT = ".exe" if IS_WINDOWS else ""
FC = environ.get("FC")
FC_REASON = "make must be used with gfortran"


def run_makefile(target):
    assert Path("makefile").is_file(), f"makefile does not exist in {os.getcwd()}"

    base_target = os.path.basename(target)
    base_message = (
        f" Rerunning {os.path.basename(__file__)} in the distribution "
        "directory will likely resolve CI failures."
    )

    # clean prior to make
    print(f"clean {base_target} with makefile")
    os.system("make clean")

    # build MODFLOW 6 with makefile
    print(f"build {base_target} with makefile")
    return_code = os.system(f"make FC={environ.get('FC', 'gfortran')}")

    assert return_code == 0, f"could not make '{base_target}'." + base_message
    assert os.path.isfile(target), f"{base_target} does not exist." + base_message


def build_mf6_makefile():
    target = "mf6"
    excludefiles = str(PROJ_ROOT_PATH / "pymake" / "excludefiles.txt")
    make_path = PROJ_ROOT_PATH / "make"
    make_path.mkdir(exist_ok=True)
    print(f"Creating makefile for {target}")
    with set_dir(make_path):
        pymake.main(
            srcdir=str(PROJ_ROOT_PATH / "src"),
            target=target,
            appdir=str(PROJ_ROOT_PATH / "bin"),
            include_subdirs=True,
            excludefiles=excludefiles,
            inplace=True,
            dryrun=True,
            makefile=True,
            networkx=True,
        )


def build_zbud6_makefile():
    target = "zbud6"
    util_path = PROJ_ROOT_PATH / "utils" / "zonebudget"
    (util_path / "make").mkdir(exist_ok=True)
    print(f"Creating makefile for {target}")
    with set_dir(util_path / "make"):
        returncode = pymake.main(
            srcdir=str(util_path / "src"),
            target=target,
            appdir=str(PROJ_ROOT_PATH / "bin"),
            extrafiles=str(util_path / "pymake" / "extrafiles.txt"),
            inplace=True,
            include_subdirs=True,
            makefile=True,
            dryrun=True,
            networkx=True,
        )

        assert returncode == 0, f"Failed to create makefile for '{target}'"


def build_mf5to6_makefile():
    target = "mf5to6"
    util_path = PROJ_ROOT_PATH / "utils" / "mf5to6"
    (util_path / "make").mkdir(exist_ok=True)
    print(f"Creating makefile for {target}")
    with set_dir(util_path / "make"):
        extrafiles = str(util_path / "pymake" / "extrafiles.txt")

        # build modflow 5 to 6 converter
        returncode = pymake.main(
            srcdir=str(util_path / "src"),
            target=target,
            appdir=str(PROJ_ROOT_PATH / "bin"),
            include_subdirs=True,
            extrafiles=extrafiles,
            inplace=True,
            dryrun=True,
            makefile=True,
            networkx=True,
            fflags=["-fall-intrinsics"],
        )

        assert returncode == 0, f"Failed to create makefile for '{target}'"


@flaky
@no_parallel
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_mf6_makefile():
    makefile_paths = [
        PROJ_ROOT_PATH / "make" / "makefile",
        PROJ_ROOT_PATH / "make" / "makedefaults",
    ]
    try:
        build_mf6_makefile()
        for p in makefile_paths:
            assert p.is_file()
    finally:
        for p in makefile_paths:
            p.unlink(missing_ok=True)


@flaky
@no_parallel
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_zbud6_makefile():
    util_path = PROJ_ROOT_PATH / "utils" / "zonebudget"
    makefile_paths = [
        util_path / "make" / "makefile",
        util_path / "make" / "makedefaults",
    ]
    try:
        build_zbud6_makefile()
        for p in makefile_paths:
            assert p.is_file()
    finally:
        for p in makefile_paths:
            p.unlink(missing_ok=True)


@flaky
@no_parallel
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_mf5to6_makefile():
    util_path = PROJ_ROOT_PATH / "utils" / "mf5to6"
    makefile_paths = [
        util_path / "make" / "makefile",
        util_path / "make" / "makedefaults",
    ]
    try:
        build_mf5to6_makefile()
        for p in makefile_paths:
            assert p.is_file()
    finally:
        for p in makefile_paths:
            p.unlink(missing_ok=True)


@flaky
@no_parallel
@requires_exe("make")
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_mf6_with_make():
    target = PROJ_ROOT_PATH / "bin" / f"mf6{EXE_EXT}"
    makefile_paths = [
        PROJ_ROOT_PATH / "make" / "makefile",
        PROJ_ROOT_PATH / "make" / "makedefaults",
    ]
    mtime = get_modified_time(target)

    try:
        build_mf6_makefile()
        with set_dir(PROJ_ROOT_PATH / "make"):
            run_makefile(target)

        # check executable was modified
        assert target.stat().st_mtime > mtime
    finally:
        print(f"clean {target} with makefile")
        os.system("make clean")
        for p in makefile_paths:
            p.unlink(missing_ok=True)


@flaky
@no_parallel
@requires_exe("make")
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_zbud6_with_make():
    target = PROJ_ROOT_PATH / "bin" / f"zbud6{EXE_EXT}"
    util_path = PROJ_ROOT_PATH / "utils" / "zonebudget"
    makefile_paths = [
        util_path / "make" / "makefile",
        util_path / "make" / "makedefaults",
    ]
    mtime = get_modified_time(target)

    try:
        build_zbud6_makefile()
        with set_dir(util_path / "make"):
            run_makefile(target)

        # check executable was modified
        assert target.stat().st_mtime > mtime
    finally:
        print(f"clean {target} with makefile")
        os.system("make clean")
        for p in makefile_paths:
            p.unlink(missing_ok=True)


@flaky
@no_parallel
@requires_exe("make")
@pytest.mark.skipif(FC == "ifort", reason=FC_REASON)
def test_build_mf5to6_with_make():
    target = PROJ_ROOT_PATH / "bin" / f"mf5to6{EXE_EXT}"
    util_path = PROJ_ROOT_PATH / "utils" / "mf5to6"
    makefile_paths = [
        util_path / "make" / "makefile",
        util_path / "make" / "makedefaults",
    ]
    mtime = get_modified_time(target)

    try:
        build_mf5to6_makefile()
        with set_dir(util_path / "make"):
            run_makefile(target)

        # check executable was modified
        assert target.stat().st_mtime > mtime
    finally:
        print(f"clean {target} with makefile")
        os.system("make clean")
        for p in makefile_paths:
            p.unlink(missing_ok=True)


if __name__ == "__main__":
    build_mf6_makefile()
    build_zbud6_makefile()
    build_mf5to6_makefile()

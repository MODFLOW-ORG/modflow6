import subprocess
import os
import argparse
import pathlib as pl
import platform
import shutil
import shlex

parser = argparse.ArgumentParser()
parser.add_argument("--compiler", type=str)
parser.add_argument("--buildtype", type=str)
parser.add_argument("--pixi", action="store_true")
parser.add_argument("action")
args = parser.parse_args()

os.environ["FC"] = args.compiler

builddir = pl.Path(f"_builddir_{platform.system()}_{args.compiler}_{args.buildtype}")
mf5to6_dir = pl.Path("utils/mf5to6")
mf5to6_builddir = mf5to6_dir / f"{builddir}"

if args.pixi:
    mf5to6_run_dir = pl.Path(".")
else:
    mf5to6_run_dir = mf5to6_dir

arg_extended = "-Dextended=false"
if os.getenv("BUILD_EXTENDED_MF6") is not None:
    if os.environ["BUILD_EXTENDED_MF6"] == '1':
        arg_extended = "-Dextended=true"

if args.action == "rebuild":
    for path in (builddir, mf5to6_builddir):
        if path.is_dir():
            shutil.rmtree(path)

if args.buildtype == "debug":
    setup_flag = ["-Dbuildtype=debug"]
else:
    setup_flag = ["-Dbuildtype=release"]

if not builddir.is_dir():
    # mf6 and zbud6
    if args.pixi:
        command = [
            "pixi",
            "run",
            "setup",
        ]
    else:
        command = [
            "meson",
            "setup",
        ]
    command += [
        str(builddir),
        "--prefix",
        os.getcwd(),
        "--libdir",
        "bin",
        arg_extended,
    ] + setup_flag
    print("Run:", shlex.join(command))
    subprocess.run(
        command,
        check=True,
    )
    # mf5to6
    if args.pixi:
        command = [
            "pixi",
            "run",
            "setup-mf5to6",
        ]
    else:
        command = [
            "meson",
            "setup",
        ]
    command += [
        str(builddir),
        "--prefix",
        os.getcwd(),
    ] + setup_flag
    print("Run:", shlex.join(command))
    subprocess.run(
        command,
        check=True,
        cwd=mf5to6_run_dir,
    )

# Remove all files from bin folder
bin_dir = os.path.join(os.getcwd(), "bin")
if os.path.isdir(bin_dir):
    for dir_entry in os.scandir(bin_dir):
        path = dir_entry.path
        if os.path.isfile(path):
            os.remove(path)

# mf6 and zbud6
if args.pixi:
    command = ["pixi", "run", "build",]
else:
    command = ["meson", "install", "-C"]    
command += [str(builddir)]
print("Run:", shlex.join(command))
subprocess.run(command, check=True)

# mf5to6
if args.pixi:
    command = ["pixi", "run", "build-mf5to6",]
else:
    command = ["meson", "install", "-C"]    
command += [str(builddir)]
print("Run:", shlex.join(command))
subprocess.run(command, check=True, cwd=mf5to6_run_dir)



import argparse
import importlib
import os
import subprocess
import sys
from pathlib import Path

import flopy
from conftest import project_root_path

dfn_path = project_root_path / "doc" / "mf6io" / "mf6ivar" / "dfn"
fpy_path = flopy.__path__[0]
print(f"flopy is installed in {fpy_path}")


def test_delete_mf6():
    pth = os.path.join(fpy_path, "mf6", "modflow")
    files = [
        entry for entry in os.listdir(pth) if os.path.isfile(os.path.join(pth, entry))
    ]
    delete_files(files, pth, exclude="mfsimulation.py")


def test_create_packages():
    # get list of files in mf6/modflow
    pth = os.path.join(fpy_path, "mf6", "modflow")
    list_files(pth)

    pth = os.path.join(fpy_path, "mf6", "utils")
    fn = "createpackages.py"

    # determine if createpackages.py exists
    fpth = os.path.join(pth, fn)
    print(f'testing if "{fpth}" exists')
    exist = os.path.isfile(fpth)
    assert exist, f'"{fpth}" does not exist'

    # run createpackages.py script
    print(f"running...{fn}")
    subprocess.check_output([sys.executable, "createpackages.py"], cwd=pth)

    # reload flopy
    print("reloading flopy")
    importlib.reload(flopy)

    # get updated list of files in mf6/modflow
    pth = os.path.join(fpy_path, "mf6", "modflow")
    list_files(pth)


def list_files(pth, exts=["py"]):
    print(f"\nLIST OF FILES IN {pth}")
    files = [
        entry for entry in os.listdir(pth) if os.path.isfile(os.path.join(pth, entry))
    ]
    idx = 0
    for fn in files:
        ext = os.path.splitext(fn)[1][1:].lower()
        if ext in exts:
            idx += 1
            print(f"    {idx:5d} - {fn}")


def delete_files(files, pth, exclude=None):
    if exclude is None:
        exclude = []
    else:
        if not isinstance(exclude, list):
            exclude = [exclude]

    for fn in files:
        if fn in exclude:
            continue
        fpth = os.path.join(pth, fn)
        try:
            print(f"removing...{fn}")
            os.remove(fpth)
        except:
            print(f"could not remove...{fn}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Update flopy from DFN files")
    parser.add_argument("-p", "--path", help="path to DFN files", default=str(dfn_path))
    args = parser.parse_args()

    path = Path(args.path).expanduser().resolve()
    print(f"Updating flopy packages from DFN files in: {path}")

    test_delete_mf6()
    test_create_packages()

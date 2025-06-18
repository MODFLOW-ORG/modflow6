import importlib

import flopy
from conftest import project_root_path
from flopy.mf6.utils import generate_classes

DFN_PATH = project_root_path / "doc" / "mf6io" / "mf6ivar" / "dfn"


if __name__ == "__main__":
    generate_classes(dfnpath=DFN_PATH, verbose=True)
    importlib.reload(flopy)

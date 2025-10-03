from pathlib import Path

from packaging.version import Version

PROJ_ROOT_PATH = Path(__file__).resolve().parent.parent
DIST_PATH = (
    PROJ_ROOT_PATH.parent
    / f"mf{Version((PROJ_ROOT_PATH / 'version.txt').read_text().strip())!s}"
)


def pytest_addoption(parser):
    # all options below are for check_dist.py
    parser.addoption("-P", "--path", action="store", default=str(DIST_PATH))
    parser.addoption("-R", "--releasemode", action="store_true", default=False)

import re
import subprocess
import sys
from collections.abc import Iterator
from datetime import datetime
from pathlib import Path
from pprint import pformat
from typing import Optional

_project_root_path = Path(__file__).resolve().parent.parent


def get_project_root_path():
    return _project_root_path


def get_modified_time(path: Path) -> float:
    return path.stat().st_mtime if path.is_file() else datetime.today().timestamp()


def glob(
    path: Path,
    pattern: str,
    included: Optional[list[str]],
    excluded: Optional[list[str]],
) -> Iterator[Path]:
    def is_included(p):
        if included is None:
            return True
        return any(i in p.name for i in included)

    def is_excluded(p):
        if excluded is None:
            return True
        return any(e in p.name for e in excluded)

    for p in path.glob(pattern):
        if p.is_file() and is_included(p) and not is_excluded(p):
            yield p


def convert_line_endings(folder, windows=True):
    """
    Convert all of the line endings to windows or unix

    """
    # Prior to zipping, enforce os line endings on all text files
    print("Converting line endings...")
    platform = sys.platform
    cmd = None
    if platform.lower() == "darwin":
        if windows:
            cmd = "find . -name '*' | xargs unix2dos"
        else:
            cmd = "find . -name '*' | xargs dos2unix"
    else:
        if windows:
            cmd = 'for /R %G in (*) do unix2dos "%G"'
        else:
            cmd = 'for /R %G in (*) do dos2unix "%G"'
    p = subprocess.Popen(cmd, cwd=folder, shell=True)
    print(p.communicate())


def split_nonnumeric(s):
    match = re.compile("[^0-9]").search(s)
    return [s[: match.start()], s[match.start() :]] if match else s


def match(l, r):
    l = set(l)
    r = set(r)
    diff = l ^ r
    return not any(diff)


def assert_match(l, r, lname=None, rname=None):
    l = set(l)
    r = set(r)
    diff = l ^ r
    lname = lname or "l"
    rname = rname or "r"
    assert not any(diff), (
        f"=> symmetric difference:\n{pformat(diff)}\n"
        f"=> {lname} - {rname}:\n{pformat(l - r)}\n"
        f"=> {rname} - {lname}:\n{pformat(r - l)}\n"
    )

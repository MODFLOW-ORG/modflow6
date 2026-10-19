"""
Apply the release-time changeset: version strings and timestamps in release
mode (update_version.py), the release's row in the release history
(../doc/ReleaseNotes/release_history.py), and formatting of the generated
version.f90. The release workflow runs this once and builds the result, and
when releasing, that same commit is what merges to master.

This is only for release builds. Development builds use the version metadata
already committed on develop.

If the release history already has a row for the version (a retried release,
or one added by hand), its date is updated and its DOI kept, so the history
and citation agree with the timestamps.
"""

import argparse
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import pytest
from packaging.version import Version

project_root_path = Path(__file__).parents[1]
distribution_path = project_root_path / "distribution"
release_notes_path = project_root_path / "doc" / "ReleaseNotes"


def prepare_release(version: str, release_date: date, doi: str | None = None):
    date_str = release_date.isoformat()

    subprocess.run(
        [
            sys.executable,
            "update_version.py",
            "-v",
            version,
            "--date",
            date_str,
            "--releasemode",
        ],
        cwd=distribution_path,
        check=True,
    )

    cmd = [
        sys.executable,
        "release_history.py",
        "-v",
        version,
        "--date",
        date_str,
        "--update-date",
    ]
    if doi:
        cmd += ["--doi", doi]
    subprocess.run(cmd, cwd=release_notes_path, check=True)

    subprocess.run(
        [
            "fprettify",
            "-c",
            str(project_root_path / ".fprettify.yaml"),
            str(project_root_path / "src" / "Utilities" / "version.f90"),
        ],
        check=True,
    )


# files the release metadata scripts read or write, and the scripts themselves,
# copied to a scratch tree so the test doesn't touch the repository
_tree_files = [
    ".fprettify.yaml",
    "version.txt",
    "meson.build",
    "pixi.toml",
    "README.md",
    "DISCLAIMER.md",
    "CITATION.cff",
    "code.json",
    ".build_rtd_docs/Doxyfile",
    "utils/mf5to6/meson.build",
    "doc/version.tex",
    "doc/version.py",
    "doc/ReleaseNotes/ReleaseNotes.tex",
    "doc/ReleaseNotes/release_history.py",
    "src/Utilities/version.f90.in",
    "src/Utilities/version.f90",
    "distribution/prepare_release.py",
    "distribution/update_version.py",
    "distribution/utils.py",
]


@pytest.fixture
def tree(tmp_path):
    for name in _tree_files:
        dst = tmp_path / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(project_root_path / name, dst)
    return tmp_path


def _snapshot(root: Path) -> dict[str, str]:
    return {name: (root / name).read_text() for name in _tree_files}


def test_prepare_release(tree):
    current = Version((tree / "version.txt").read_text().strip())
    version = f"{current.major}.{current.minor + 1}.0"
    doi = "https://doi.org/10.5066/TEST"
    cmd = [
        sys.executable,
        str(tree / "distribution" / "prepare_release.py"),
        "-v",
        version,
        "--date",
        "2026-12-01",
        "--doi",
        doi,
    ]

    subprocess.run(cmd, cwd=tree / "distribution", check=True)
    first = _snapshot(tree)

    assert first["version.txt"].strip() == version
    assert "IDEVELOPMODE = 0" in first["src/Utilities/version.f90"]
    assert (
        f"{version} & December 1, 2026 & \\url{{{doi}}}"
        in first["doc/ReleaseNotes/ReleaseNotes.tex"]
    )

    # rerunning with the same inputs changes nothing but the run time embedded
    # in doc/version.py, which the release workflow relies on to skip the commit
    subprocess.run(cmd, cwd=tree / "distribution", check=True)
    second = _snapshot(tree)
    del first["doc/version.py"], second["doc/version.py"]
    assert first == second


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        required=True,
        help="Version number of the release, e.g. 6.9.0",
    )
    parser.add_argument(
        "--date",
        type=date.fromisoformat,
        default=None,
        help="Date of the release (YYYY-MM-DD). Defaults to today.",
    )
    parser.add_argument(
        "--doi",
        default=None,
        help="DOI of the release, bare or as a link. Needed for a new minor "
        "release unless its release history row already exists. Defaults to "
        "the DOI of the latest release with the same major and minor version.",
    )
    args = parser.parse_args()

    try:
        prepare_release(args.version, args.date or date.today(), doi=args.doi)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

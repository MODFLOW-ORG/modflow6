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
import subprocess
import sys
from datetime import date
from pathlib import Path

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

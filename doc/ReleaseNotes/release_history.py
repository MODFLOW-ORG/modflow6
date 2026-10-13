"""Add a release to the release history table in ReleaseNotes.tex.

The table lists the version, date and DOI of each release, e.g.

    6.8.0 & September 2, 2026 & \\url{https://doi.org/10.5066/P1PGE9XW} \\\\

Adding a release appends a row for the version, unless the table already has one,
in which case it is left as is, or with --update-date only its date is updated. The
release workflow does this, so a release retried on a later day, which is a new
build, is dated by that build everywhere. The date is today's unless --date is
given. The DOI is --doi if given, otherwise that of the latest release with the
same major and minor version: DOIs change with minor releases, not patch
releases, so a patch release reuses the DOI of the minor release it patches. A new
minor release needs --doi, unless its row was already added by hand.

Only release versions (major.minor.patch) are added. Release candidates and
development versions are skipped.

Use --check to verify a row could be added without changing anything. This script
only needs the standard library, so it can run without an environment.

The version is read from version.txt if --version isn't given. The table is the
source of a release's date and DOI: update_version.py reads them from it for the
software citation, and mk_releasenotes.py and reset_releasenotes.py for the archive
header of a release's notes.
"""

import argparse
import re
import sys
from datetime import date, datetime
from pathlib import Path

notes_dir = Path(__file__).parent
notes_path = notes_dir / "ReleaseNotes.tex"
version_path = notes_dir.parents[1] / "version.txt"

# a row of the table, e.g.
# 6.8.0 & September 2, 2026 & \url{https://doi.org/10.5066/P1PGE9XW} \\
_row = re.compile(
    r"^[ \t]*(?P<version>\d+\.\d+\.\d+)[ \t]*&[ \t]*(?P<date>[^&\n]+?)[ \t]*&"
    r"[ \t]*\\url\{(?P<doi>[^}\n]+)\}[ \t]*\\\\[ \t]*$",
    re.MULTILINE,
)

# a DOI, bare or as a link. Only characters found in DOIs and safe in LaTeX.
_doi = re.compile(
    r"(?:https?://(?:dx\.)?doi\.org/)?(?P<doi>10\.\d{4,9}/[A-Za-z0-9._;()/:-]+)"
)


def is_release(version: str) -> bool:
    """Whether the version is a release version: major.minor.patch, no suffix."""
    return re.fullmatch(r"\d+\.\d+\.\d+", version) is not None


def normalize_doi(doi: str) -> str:
    """Canonical link for a DOI given bare (10.5066/P1PGE9XW) or as a link.
    Raises a ValueError if it isn't recognizable as a DOI."""
    match = _doi.fullmatch(doi.strip())
    if not match:
        raise ValueError(
            f"Invalid DOI {doi!r}, expected e.g. https://doi.org/10.5066/P1PGE9XW"
        )
    return f"https://doi.org/{match['doi']}"


def _rows(path: Path) -> list[re.Match]:
    rows = list(_row.finditer(path.read_text()))
    if not rows:
        raise ValueError(f"No release history table rows found in {path}")
    return rows


def latest_release(path: Path = notes_path) -> tuple[str, str]:
    """Version and date of the most recent release, from the last table row."""
    row = _rows(path)[-1]
    return row["version"], row["date"]


def find_release(version: str, path: Path = notes_path) -> tuple[date, str] | None:
    """Date and DOI in the table row for the version, or None if it has no row.
    A ValueError is raised if the row's date isn't in the format 'April 8, 2026'."""
    for row in _rows(path):
        if row["version"] == version:
            try:
                release_date = datetime.strptime(row["date"], "%B %d, %Y").date()
            except ValueError:
                raise ValueError(
                    f"Unrecognized date {row['date']!r} in the row for {version} "
                    f"in {path}, expected e.g. 'April 8, 2026'"
                ) from None
            return release_date, row["doi"]
    return None


def _format_date(release_date: date) -> str:
    return f"{release_date:%B} {release_date.day}, {release_date.year}"


def add_release(
    version: str,
    release_date: date,
    doi: str | None = None,
    path: Path = notes_path,
    dry_run: bool = False,
) -> bool:
    """
    Add a row for the version after the last row of the table. Returns True if a
    row was added, False if the table already has a row for the version.

    If no DOI is given, the DOI of the latest release with the same major and minor
    version is used. A ValueError is raised if there is none, if the DOI is
    invalid, if the version isn't a release version, or if there is no table.

    With dry_run, nothing is written, but the same errors are raised.
    """

    if not is_release(version):
        raise ValueError(f"{version} is not a release version")

    rows = _rows(path)
    if any(row["version"] == version for row in rows):
        return False

    if doi:
        doi = normalize_doi(doi)
    else:
        series = version.rsplit(".", 1)[0]
        same_series = [r for r in rows if r["version"].rsplit(".", 1)[0] == series]
        if not same_series:
            raise ValueError(
                f"No DOI for new release {version}, and no previous {series}.x "
                f"release to take one from. Provide a DOI, or add a row for "
                f"{version} to the release history in {path.name}."
            )
        doi = same_series[-1]["doi"]

    if not dry_run:
        text = path.read_text()
        row = f"{version} & {_format_date(release_date)} & \\url{{{doi}}} \\\\"
        end = rows[-1].end()
        path.write_text(text[:end] + "\n" + row + text[end:])
        print(f"Added release history row to {path}: {row}", file=sys.stderr)
    return True


def update_release_date(
    version: str,
    release_date: date,
    path: Path = notes_path,
    dry_run: bool = False,
) -> bool:
    """
    Set the date in the table row for the version, keeping its DOI. Returns True if
    the date changed, False if it already had that date. A ValueError is raised if
    the table has no row for the version.

    With dry_run, nothing is written, but the same errors are raised.
    """

    row = next((r for r in _rows(path) if r["version"] == version), None)
    if row is None:
        raise ValueError(f"No release history row for {version} in {path}")

    date_str = _format_date(release_date)
    if row["date"] == date_str:
        return False

    if not dry_run:
        text = path.read_text()
        start, end = row.span("date")
        path.write_text(text[:start] + date_str + text[end:])
        print(
            f"Updated date of release history row for {version} in {path}: "
            f"{row['date']} -> {date_str}",
            file=sys.stderr,
        )
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        default=None,
        help="Version of the release. Defaults to the contents of version.txt.",
    )
    parser.add_argument(
        "--date",
        type=date.fromisoformat,
        default=None,
        help="Release date (YYYY-MM-DD). Defaults to today.",
    )
    parser.add_argument(
        "-d",
        "--doi",
        default=None,
        help="DOI of the release, bare or as a link, e.g. "
        "https://doi.org/10.5066/P1PGE9XW. Defaults to the DOI of the latest "
        "release with the same major and minor version.",
    )
    parser.add_argument(
        "--update-date",
        default=False,
        action="store_true",
        help="If the table already has a row for the version, update its date "
        "instead of leaving it as is. Its DOI is kept.",
    )
    parser.add_argument(
        "--check",
        default=False,
        action="store_true",
        help="Check that a row could be added, without changing anything. "
        "Exits with an error if not.",
    )
    parser.add_argument(
        "--notes",
        type=Path,
        default=notes_path,
        help="Release notes file with the table. Defaults to ReleaseNotes.tex.",
    )
    args = parser.parse_args()

    version = args.version or version_path.read_text().strip()
    if not is_release(version):
        print(f"{version} is not a release version, skipping", file=sys.stderr)
        sys.exit(0)

    release_date = args.date or date.today()
    try:
        added = add_release(
            version,
            release_date,
            doi=args.doi,
            path=args.notes,
            dry_run=args.check,
        )
        updated = (
            not added
            and args.update_date
            and update_release_date(
                version, release_date, path=args.notes, dry_run=args.check
            )
        )
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(1)

    if not added and not updated:
        print(f"{args.notes} already has a row for {version}", file=sys.stderr)
    elif args.check:
        print(f"A row for {version} can be added to {args.notes}", file=sys.stderr)

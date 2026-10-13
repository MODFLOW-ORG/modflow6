import subprocess
import sys
from datetime import date
from pathlib import Path

import pytest
from release_history import (
    _row,
    add_release,
    find_release,
    is_release,
    latest_release,
    normalize_doi,
    notes_path,
    update_release_date,
)

SCRIPT = Path(__file__).parent / "release_history.py"

SAMPLE = r"""\begin{tabular*}{\columnwidth}{l l l}
6.7.0 & February 6, 2026 & \url{https://doi.org/10.5066/P1IJAXDZ} \\
6.8.0 & September 2, 2026 & \url{https://doi.org/10.5066/P1PGE9XW} \\
\hline
\label{tab:releases}
\end{tabular*}
"""


@pytest.fixture
def notes(tmp_path):
    path = tmp_path / "ReleaseNotes.tex"
    path.write_text(SAMPLE)
    return path


def test_add_release_patch_reuses_doi(notes):
    assert add_release("6.8.1", date(2026, 9, 21), path=notes)

    expected = SAMPLE.replace(
        "\\hline",
        "6.8.1 & September 21, 2026 & \\url{https://doi.org/10.5066/P1PGE9XW} \\\\\n"
        "\\hline",
    )
    assert notes.read_text() == expected


@pytest.mark.parametrize(
    "doi",
    ["https://doi.org/10.5066/NEW1", "10.5066/NEW1", " http://doi.org/10.5066/NEW1 "],
)
def test_add_release_minor_with_doi(notes, doi):
    assert add_release("6.9.0", date(2026, 12, 1), doi=doi, path=notes)

    lines = notes.read_text().splitlines()
    assert (
        lines[3] == r"6.9.0 & December 1, 2026 & \url{https://doi.org/10.5066/NEW1} \\"
    )
    assert lines[4] == "\\hline"


def test_add_release_minor_without_doi(notes):
    with pytest.raises(ValueError, match=r"No DOI for new release 6\.9\.0"):
        add_release("6.9.0", date(2026, 12, 1), path=notes)
    assert notes.read_text() == SAMPLE


def test_add_release_series_is_exact(notes):
    # 6.8.x must not be mistaken for 6.80.x or 6.1.x for 6.10.x
    with pytest.raises(ValueError, match="No DOI"):
        add_release("6.80.0", date(2026, 12, 1), path=notes)
    with pytest.raises(ValueError, match="No DOI"):
        add_release("6.70.1", date(2026, 12, 1), path=notes)


def test_add_release_existing_row_is_kept(notes):
    assert add_release("6.8.1", date(2026, 9, 21), path=notes)
    added = notes.read_text()

    # same version again, and a version from before
    assert not add_release("6.8.1", date(2026, 9, 22), path=notes)
    assert not add_release("6.8.0", date(2026, 9, 22), doi="10.5066/OTHER", path=notes)
    assert notes.read_text() == added


def test_add_release_dry_run(notes):
    assert add_release("6.8.1", date(2026, 9, 21), path=notes, dry_run=True)
    assert notes.read_text() == SAMPLE

    with pytest.raises(ValueError, match="No DOI"):
        add_release("6.9.0", date(2026, 12, 1), path=notes, dry_run=True)


@pytest.mark.parametrize(
    "version", ["6.9.0rc", "6.9.0rc1", "6.9.0.dev0", "6.9", "v6.9.0"]
)
def test_add_release_not_a_release_version(notes, version):
    assert not is_release(version)
    with pytest.raises(ValueError, match="not a release version"):
        add_release(version, date(2026, 12, 1), doi="10.5066/NEW1", path=notes)


@pytest.mark.parametrize(
    "doi",
    [
        "",
        "P1PGE9XW",
        "https://example.com/10.5066/P1PGE9XW",
        "https://doi.org/10.5066/A}B",
        "https://doi.org/10.5066/A B",
        "https://doi.org/10.5066/A\\B",
        "https://doi.org/10.5066/A%B",
    ],
)
def test_normalize_doi_invalid(doi):
    with pytest.raises(ValueError, match="Invalid DOI"):
        normalize_doi(doi)


def test_add_release_invalid_doi(notes):
    with pytest.raises(ValueError, match="Invalid DOI"):
        add_release("6.9.0", date(2026, 12, 1), doi="oops}", path=notes)
    assert notes.read_text() == SAMPLE


def test_update_release_date(notes):
    assert add_release("6.8.1", date(2026, 9, 21), path=notes)
    added = notes.read_text()

    # same date: nothing to do
    assert not update_release_date("6.8.1", date(2026, 9, 21), path=notes)
    assert notes.read_text() == added

    # dry run: reports the change without making it
    assert update_release_date("6.8.1", date(2026, 9, 23), path=notes, dry_run=True)
    assert notes.read_text() == added

    # new date: only that row's date changes, the DOI is kept
    assert update_release_date("6.8.1", date(2026, 9, 23), path=notes)
    assert find_release("6.8.1", notes) == (
        date(2026, 9, 23),
        "https://doi.org/10.5066/P1PGE9XW",
    )
    assert notes.read_text() == added.replace("September 21", "September 23")


def test_update_release_date_no_row(notes):
    with pytest.raises(ValueError, match=r"No release history row for 6\.8\.1"):
        update_release_date("6.8.1", date(2026, 9, 21), path=notes)


@pytest.mark.parametrize("update_date", [False, True])
def test_cli_existing_row(notes, update_date):
    def run(day):
        args = [sys.executable, SCRIPT, "-v", "6.8.1", "--date", f"2026-09-{day}"]
        args += ["--notes", notes]
        if update_date:
            args.append("--update-date")
        subprocess.run(args, check=True)

    run(21)
    run(23)
    expected = date(2026, 9, 23 if update_date else 21)
    assert find_release("6.8.1", notes)[0] == expected


def test_find_release(notes):
    assert find_release("6.7.0", notes) == (
        date(2026, 2, 6),
        "https://doi.org/10.5066/P1IJAXDZ",
    )
    assert find_release("6.8.0", notes) == (
        date(2026, 9, 2),
        "https://doi.org/10.5066/P1PGE9XW",
    )
    assert find_release("6.8.1", notes) is None
    assert find_release("6.9.0.dev0", notes) is None

    # the row added for a patch release is found, with the date and reused DOI
    add_release("6.8.1", date(2026, 9, 21), path=notes)
    assert find_release("6.8.1", notes) == (
        date(2026, 9, 21),
        "https://doi.org/10.5066/P1PGE9XW",
    )


def test_find_release_bad_date(tmp_path):
    path = tmp_path / "ReleaseNotes.tex"
    path.write_text(
        "6.8.0 & 2026-09-02 & \\url{https://doi.org/10.5066/P1PGE9XW} \\\\\n"
    )
    with pytest.raises(ValueError, match="Unrecognized date"):
        find_release("6.8.0", path)


def test_latest_release(notes):
    assert latest_release(notes) == ("6.8.0", "September 2, 2026")
    add_release("6.8.1", date(2026, 9, 21), path=notes)
    assert latest_release(notes) == ("6.8.1", "September 21, 2026")


def test_no_table(tmp_path):
    path = tmp_path / "ReleaseNotes.tex"
    path.write_text("nothing to see here\n")
    with pytest.raises(ValueError, match="No release history table rows"):
        add_release("6.8.1", date(2026, 9, 21), path=path)
    with pytest.raises(ValueError, match="No release history table rows"):
        latest_release(path)


def test_repo_release_notes():
    # the table in the repository's release notes must be parseable, with
    # dates in the format the rest of the release tooling expects
    rows = list(_row.finditer(notes_path.read_text()))
    assert rows
    assert latest_release()[0] == rows[-1]["version"]
    for row in rows:
        assert find_release(row["version"]) is not None

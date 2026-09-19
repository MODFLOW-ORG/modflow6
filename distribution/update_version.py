#!/usr/bin/python

"""
Update files in this modflow6 repository according to release information.

This script is used to update several files in the modflow6 repository, including:

  ../version.txt
  ../meson.build
  ../utils/mf5to6/meson.build
  ../doc/version.tex
  ../README.md
  ../DISCLAIMER.md
  ../code.json
  ../src/Utilities/version.f90.in
  ../src/Utilities/version.f90
  ../doc/ReleaseNotes/ReleaseNotes.tex

Information in these files include version number (major.minor.patch[label]), build
timestamp, whether or not the release is preliminary/provisional or official/approved,
whether the source code should be compiled in develop mode (IDEVELOPMODE = 1) or for
release, and other metadata.

The version number is read from ../version.txt, which contains major, minor, and patch
version numbers, and an optional label. Version numbers are substituted into source
code, latex files, markdown files, etc. The version number can be provided explicitly
using --version, short -v.

If --releasemode is provided, IDEVELOPMODE is set to 0 in src/Utilities/version.f90.in
and src/Utilities/version.f90. Otherwise, IDEVELOPMODE is set to 1.

if --releasemode is provided, the disclaimer in src/Utilities/version.f90.in and the
README/DISCLAIMER markdown files is modified to reflect review and approval.
Otherwise the language reflects preliminary/provisional status.

For a release version (i.e. not a pre-release or development version), a row is added
to the release history table in ReleaseNotes.tex if there isn't one for the version
already. The row's date is today's date unless --date is provided. The DOI is taken
from --doi if provided, otherwise from the latest release with the same major and
minor version (patch releases keep the DOI of the minor release they patch). The DOI
of a new minor or major release must be provided, either via --doi or by adding the
row manually before running this script.
"""

import argparse
import json
import os
import re
import sys
import textwrap
from collections import OrderedDict
from datetime import date, datetime
from pathlib import Path

import pytest
import yaml
from filelock import FileLock
from modflow_devtools.markers import no_parallel
from packaging.version import Version

from utils import get_modified_time

project_name = "MODFLOW 6"
project_root_path = Path(__file__).resolve().parent.parent
version_file_path = project_root_path / "version.txt"
touched_file_paths = [
    version_file_path,
    project_root_path / "meson.build",
    project_root_path / "utils" / "mf5to6" / "meson.build",
    project_root_path / "doc" / "version.tex",
    project_root_path / "doc" / "version.py",
    project_root_path / "README.md",
    project_root_path / "DISCLAIMER.md",
    project_root_path / "CITATION.cff",
    project_root_path / "code.json",
    project_root_path / "src" / "Utilities" / "version.f90.in",
    project_root_path / "src" / "Utilities" / "version.f90",
]
release_notes_path = project_root_path / "doc" / "ReleaseNotes" / "ReleaseNotes.tex"

# row in the release history table in release notes, e.g.
# 6.8.0 & September 2, 2026 & \url{https://doi.org/10.5066/P1PGE9XW} \\
_release_history_row = re.compile(
    r"^[ \t]*(?P<version>\d+\.\d+\.\d+)[ \t]*&[ \t]*(?P<date>[^&\n]+?)[ \t]*&"
    r"[ \t]*\\url\{(?P<doi>[^}\n]+)\}[ \t]*\\\\[ \t]*$",
    re.MULTILINE,
)


_approved_fmtdisclaimer = '''  character(len=*), parameter :: FMTDISCLAIMER = &
    "(/,&
    &'This software has been approved for release by the U.S. Geological ',/,&
    &'Survey (USGS). Although the software has been subjected to rigorous ',/,&
    &'review, the USGS reserves the right to update the software as needed ',/,&
    &'pursuant to further analysis and review. No warranty, expressed or ',/,&
    &'implied, is made by the USGS or the U.S. Government as to the ',/,&
    &'functionality of the software and related material nor shall the ',/,&
    &'fact of release constitute any such warranty. Furthermore, the ',/,&
    &'software is released on condition that neither the USGS nor the U.S. ',/,&
    &'Government shall be held liable for any damages resulting from its ',/,&
    &'authorized or unauthorized use. Also refer to the USGS Water ',/,&
    &'Resources Software User Rights Notice for complete use, copyright, ',/,&
    &'and distribution information.',/)"'''

_preliminary_fmtdisclaimer = '''  character(len=*), parameter :: FMTDISCLAIMER = &
    "(/,&
    &'This software is preliminary or provisional and is subject to ',/,&
    &'revision. It is being provided to meet the need for timely best ',/,&
    &'science. The software has not received final approval by the U.S. ',/,&
    &'Geological Survey (USGS). No warranty, expressed or implied, is made ',/,&
    &'by the USGS or the U.S. Government as to the functionality of the ',/,&
    &'software and related material nor shall the fact of release ',/,&
    &'constitute any such warranty. The software is provided on the ',/,&
    &'condition that neither the USGS nor the U.S. Government shall be held ',/,&
    &'liable for any damages resulting from the authorized or unauthorized ',/,&
    &'use of the software.',/)"'''

_approved_disclaimer = """Disclaimer
----------

This software has been approved for release by the U.S. Geological Survey
(USGS). Although the software has been subjected to rigorous review, the USGS
reserves the right to update the software as needed pursuant to further analysis
and review. No warranty, expressed or implied, is made by the USGS or the U.S.
Government as to the functionality of the software and related material nor
shall the fact of release constitute any such warranty. Furthermore, the
software is released on condition that neither the USGS nor the U.S. Government
shall be held liable for any damages resulting from its authorized or
unauthorized use.
"""

_preliminary_disclaimer = """Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.
"""


def get_disclaimer(developmode: bool = False, formatted: bool = False) -> str:
    if developmode:
        return _preliminary_fmtdisclaimer if formatted else _preliminary_disclaimer
    return _approved_fmtdisclaimer if formatted else _approved_disclaimer


# Umbrella DOI for the "MODFLOW and Related Programs" software release page.
# Used as the software citation DOI when a release-specific one isn't passed
# via --doi or found in the release history.
_default_doi = "https://doi.org/10.5066/F76Q1VQV"


def get_software_citation(
    timestamp: datetime,
    version: Version,
    doi: str = _default_doi,
    developmode: bool = False,
) -> str:
    # get data Software/Code citation for FloPy
    citation = yaml.safe_load((project_root_path / "CITATION.cff").read_text())

    # format author names
    authors = []
    for author in citation["authors"]:
        tauthor = author["family-names"] + ", "
        gnames = author["given-names"].split()
        if len(gnames) > 1:
            for gname in gnames:
                tauthor += gname[0]
                if len(gname) > 1:
                    tauthor += "."
                # tauthor += " "
        else:
            tauthor += author["given-names"]
        authors.append(tauthor.rstrip())

    line = ""
    for ipos, tauthor in enumerate(authors):
        if ipos > 0:
            line += ", "
        if ipos == len(authors) - 1:
            line += "and "
        # add formatted author name to line
        line += tauthor

    # add the rest of the citation
    line += (
        f", {timestamp.year}, "
        f"MODFLOW 6 Modular Hydrologic Model version {version}: "
        f"U.S. Geological Survey Software Release, {timestamp:%-d %B %Y}, "
        f"{doi}"
    )

    return line


def log_update(path, version: Version):
    print(f"Updated {path} with version {version}", file=sys.stderr)


def update_version_txt_and_py(version: Version, timestamp: datetime):
    with open(version_file_path, "w") as f:
        f.write(str(version))
    log_update(version_file_path, version)

    py_path = project_root_path / "doc" / version_file_path.name.replace(".txt", ".py")
    with open(py_path, "w") as f:
        f.write(
            f"# {project_name} version file automatically "
            + f"created using...{os.path.basename(__file__)}\n"
        )
        f.write("# created on..." + f"{timestamp.strftime('%B %d, %Y %H:%M:%S')}\n")
        f.write(f'__version__ = "{version}"\n')
    log_update(py_path, version)


def update_meson_build(version: Version):
    paths = [
        project_root_path / "meson.build",
        project_root_path / "utils" / "mf5to6" / "meson.build",
    ]
    for path in paths:
        lines = open(path, "r").read().splitlines()
        with open(path, "w") as f:
            for line in lines:
                if "version:" in line and "meson_version:" not in line:
                    line = f"  version: '{version}',"
                f.write(f"{line}\n")
        log_update(path, version)


def update_version_tex(version: Version, timestamp: datetime, developmode: bool = True):
    path = project_root_path / "doc" / "version.tex"
    with open(path, "w") as f:
        lines = [
            "\\newcommand{\\modflowversion}{mf" + str(version) + "}",
            "\\newcommand{\\modflowdate}{" + f"{timestamp.strftime('%B %d, %Y')}" + "}",
            (
                "\\newcommand{\\currentmodflowversion} "
                "{Version \\modflowversion---\\modflowdate}"
            ),
            "\\newif\\ifdevelopmode",
            f"\\developmode{'true' if developmode else 'false'}",
        ]
        for line in lines:
            f.write(f"{line}\n")

    log_update(path, version)


def update_version_f90(
    version: Version | None,
    timestamp: datetime,
    developmode: bool = False,
):
    version_spl = str(version).rpartition("-")
    version_num = version_spl[0] if version_spl[1] else str(version_spl[2])
    new_title = "" if developmode else f" {timestamp.strftime('%m/%d/%Y')}"

    template_path = project_root_path / "src" / "Utilities" / "version.f90.in"
    static_path = project_root_path / "src" / "Utilities" / "version.f90"

    lines = open(template_path, "r").read().splitlines()
    updated_lines = []
    skip = False
    for line in lines:
        if skip:
            if ',/)"' in line:
                skip = False
            continue
        elif ":: IDEVELOPMODE =" in line:
            line = (
                "  integer(I4B), parameter :: "
                + f"IDEVELOPMODE = {1 if developmode else 0}"
            )
        elif ":: VERSIONNUMBER =" in line:
            line = line.rpartition("::")[0] + f":: VERSIONNUMBER = '{version_num}'"
        elif ":: VERSIONVCSTAG =" in line and not developmode:
            # release builds run before the release commit/tag exists,
            # so set an empty tag here rather than rely on meson at build time
            line = line.replace("@VCS_TAG@", "")
        elif ":: VERSIONTITLE =" in line:
            line = line.rpartition("::")[0] + f":: VERSIONTITLE = '{new_title}'"
        elif ":: FMTDISCLAIMER =" in line:
            line = get_disclaimer(developmode=developmode, formatted=True)
            skip = True
        updated_lines.append(line)

    with open(template_path, "w") as f:
        for line in updated_lines:
            f.write(f"{line}\n")
    log_update(template_path, version)

    with open(static_path, "w") as f:
        for line in updated_lines:
            f.write(f"{line.replace('@VCS_TAG@', '')}\n")
    log_update(static_path, version)


def update_readme_and_disclaimer(version: Version, developmode: bool = False):
    disclaimer = get_disclaimer(developmode, formatted=False)
    readme_path = str(project_root_path / "README.md")
    readme_lines = open(readme_path, "r").read().splitlines()
    with open(readme_path, "w") as f:
        for line in readme_lines:
            if "## Version " in line:
                f.write(f"### Version {version}\n")
            elif "Disclaimer" in line:
                f.write(f"{disclaimer}\n")
                break
            else:
                f.write(f"{line}\n")
    log_update(readme_path, version)

    disclaimer_path = project_root_path / "DISCLAIMER.md"
    with open(disclaimer_path, "w") as f:
        f.write(disclaimer)
    log_update(disclaimer_path, version)


def update_citation_cff(version: Version, timestamp: datetime):
    path = project_root_path / "CITATION.cff"
    citation = yaml.safe_load(path.read_text())
    citation["version"] = str(version)
    citation["date-released"] = timestamp.strftime("%Y-%m-%d")

    with open(path, "w") as f:
        yaml.safe_dump(
            citation, f, allow_unicode=True, default_flow_style=False, sort_keys=False
        )
    log_update(path, version)


def update_codejson(version: Version, timestamp: datetime, developmode: bool = False):
    path = project_root_path / "code.json"
    with open(path, "r") as f:
        data = json.load(f, object_pairs_hook=OrderedDict)

    data[0]["date"]["metadataLastUpdated"] = timestamp.strftime("%Y-%m-%d")
    data[0]["version"] = str(version)
    data[0]["status"] = "Preliminary" if developmode else "Release"
    with open(path, "w") as f:
        json.dump(data, f, indent=4)
        f.write("\n")

    log_update(path, version)


def update_doxyfile(version: Version):
    path = project_root_path / ".build_rtd_docs" / "Doxyfile"
    lines = open(path, "r").readlines()
    tag = "PROJECT_NUMBER"
    with open(path, "w") as fp:
        for line in lines:
            if tag in line:
                line = f'{tag}         = "version {version}"\n'
            fp.write(line)


def update_pixi(version: Version):
    path = project_root_path / "pixi.toml"
    lines = open(path, "r").readlines()
    tag = "version ="
    with open(path, "w") as fp:
        for line in lines:
            if line.startswith(tag):
                line = f'{tag} "{version}"\n'
            fp.write(line)


def get_release_history_doi(
    version: Version, path: Path = release_notes_path
) -> str | None:
    """
    DOI in the release history table in the release notes for the given
    version, or None if the table has no row for the version.
    """

    ver = f"{version.major}.{version.minor}.{version.micro}"
    for row in _release_history_row.finditer(path.read_text()):
        if row["version"] == ver:
            return row["doi"]
    return None


def update_release_history(
    version: Version,
    release_date: date,
    doi: str | None = None,
    path: Path = release_notes_path,
) -> bool:
    """
    Add a row for the given version to the release history table in the release
    notes, after the last existing row. Returns True if a row was added, False
    if the table already has a row for the version, which is left as is.

    If no DOI is provided, the DOI of the latest release with the same major and
    minor version is reused, as DOIs change with minor releases but not patches.
    A ValueError is raised if there is no such release, or no table.
    """

    text = path.read_text()
    rows = list(_release_history_row.finditer(text))
    if not rows:
        raise ValueError(f"No release history table rows found in {path}")

    ver = f"{version.major}.{version.minor}.{version.micro}"
    if any(row["version"] == ver for row in rows):
        print(f"{path} already has a release history row for {ver}", file=sys.stderr)
        return False

    if doi is None:
        series = [
            row
            for row in rows
            if row["version"].split(".")[:2] == [str(version.major), str(version.minor)]
        ]
        if not series:
            raise ValueError(
                f"No DOI for new release {ver}, and no previous {version.major}."
                f"{version.minor}.x release to take one from. Pass a DOI with "
                f"--doi, or add a row for {ver} to the release history in {path}."
            )
        doi = series[-1]["doi"]

    date_str = f"{release_date:%B} {release_date.day}, {release_date.year}"
    row = f"{ver} & {date_str} & \\url{{{doi}}} \\\\"
    text = text[: rows[-1].end()] + "\n" + row + text[rows[-1].end() :]
    path.write_text(text)
    print(f"Added release history row to {path}: {row}", file=sys.stderr)
    return True


def update_version(
    version: Version = None,
    timestamp: datetime = datetime.now(),
    developmode: bool = False,
    release_date: date | None = None,
    doi: str | None = None,
):
    """
    Update version information stored in version.txt in the project root,
    as well as several other files in the repository. Version updates are
    performed by explicitly providing a version argument to this function
    and a lock is held on the version file to make sure that the state of
    the multiple files containing version information stays synchronized.
    If no version argument is provided, the version number isn't changed.

    If the version is a release version (not a pre-release or development
    version), a row is added to the release history in the release notes, if
    there isn't one already, dated with release_date (default: the timestamp's
    date) and linking to doi (default: the DOI of the minor release, see
    update_release_history).
    """

    lock_path = Path(version_file_path.name + ".lock")
    try:
        lock = FileLock(lock_path)
        previous = Version(version_file_path.read_text().strip())
        version = version if version else previous

        with lock:
            # first, since it can fail and nothing should be modified if it does
            if not (version.is_prerelease or version.is_devrelease):
                update_release_history(version, release_date or timestamp.date(), doi)
            update_version_txt_and_py(version, timestamp)
            update_meson_build(version)
            update_version_tex(version, timestamp, developmode)
            update_version_f90(version, timestamp, developmode)
            update_readme_and_disclaimer(version, developmode)
            update_citation_cff(version, timestamp)
            update_codejson(version, timestamp, developmode)
            update_doxyfile(version)
            update_pixi(version)

    finally:
        lock_path.unlink(missing_ok=True)


_initial_version = Version("0.0.1")
_current_version = Version(version_file_path.read_text().strip())


def release_version() -> Version:
    """Current development version, any development segment (e.g. '.dev0') removed."""
    return Version(_current_version.base_version)


def post_release_version() -> Version:
    """Development version for the next cycle: minor incremented, '.dev0' suffix."""
    version = Version(_current_version.base_version)
    return Version(f"{version.major}.{version.minor + 1}.0.dev0")


@no_parallel
@pytest.mark.skip(reason="reverts repo files on cleanup, treat carefully")
@pytest.mark.parametrize(
    "version",
    [
        None,
        _initial_version,
        Version(
            f"{_initial_version.major}.{_initial_version.minor}.dev{_initial_version.micro}"
        ),
    ],
)
@pytest.mark.parametrize("full", [True, False])
def test_update_version(version, full):
    m_times = [get_modified_time(file) for file in touched_file_paths]
    timestamp = datetime.now()

    try:
        update_version(
            timestamp=timestamp,
            version=version,
            developmode=full,
            doi="https://doi.org/10.5066/TEST",
        )
        updated = Version(version_file_path.read_text().strip())

        # check files containing version info were modified
        for p, t in zip(touched_file_paths, m_times):
            assert p.stat().st_mtime > t

        # check version number and optional label are correct
        if version:
            # version should be auto-incremented
            assert updated == _initial_version
        else:
            # version should not have changed
            assert updated == _current_version

        # check IDEVELOPMODE was set correctly
        version_f90_path = project_root_path / "src" / "Utilities" / "version.f90.in"
        lines = version_f90_path.read_text().splitlines()
        assert any(f"IDEVELOPMODE = {0 if full else 1}" in line for line in lines)

        # check disclaimer has appropriate language
        disclaimer_path = project_root_path / "DISCLAIMER.md"
        lines = disclaimer_path.read_text().splitlines()
        assert any(("approved for release") in line for line in lines) == full
        assert any(("preliminary or provisional") in line for line in lines) != full

    finally:
        for p in [*touched_file_paths, release_notes_path]:
            os.system(f"git restore {p}")


_release_history_sample = r"""\begin{tabular*}{\columnwidth}{l l l}
6.7.0 & February 6, 2026 & \url{https://doi.org/10.5066/P1IJAXDZ} \\
6.8.0 & September 2, 2026 & \url{https://doi.org/10.5066/P1PGE9XW} \\
\hline
\label{tab:releases}
\end{tabular*}
"""


@pytest.fixture
def release_notes(tmp_path):
    path = tmp_path / "ReleaseNotes.tex"
    path.write_text(_release_history_sample)
    return path


def test_update_release_history_patch_reuses_doi(release_notes):
    assert update_release_history(
        Version("6.8.1"), date(2026, 9, 21), path=release_notes
    )

    expected = _release_history_sample.replace(
        "\\hline",
        "6.8.1 & September 21, 2026 & \\url{https://doi.org/10.5066/P1PGE9XW} \\\\\n"
        "\\hline",
    )
    assert release_notes.read_text() == expected


def test_update_release_history_minor_with_doi(release_notes):
    doi = "https://doi.org/10.5066/NEWDOI"
    assert update_release_history(
        Version("6.9.0"), date(2026, 12, 1), doi=doi, path=release_notes
    )

    lines = release_notes.read_text().splitlines()
    assert lines[3] == f"6.9.0 & December 1, 2026 & \\url{{{doi}}} \\\\"
    assert lines[4] == "\\hline"


def test_update_release_history_minor_without_doi(release_notes):
    with pytest.raises(ValueError, match=r"No DOI for new release 6\.9\.0"):
        update_release_history(Version("6.9.0"), date(2026, 12, 1), path=release_notes)
    assert release_notes.read_text() == _release_history_sample


def test_update_release_history_existing_row_is_kept(release_notes):
    # e.g. a minor release whose row was added manually, or a repeated run
    assert update_release_history(
        Version("6.8.1"), date(2026, 9, 21), path=release_notes
    )
    added = release_notes.read_text()
    assert not update_release_history(
        Version("6.8.1"), date(2026, 9, 22), path=release_notes
    )
    assert release_notes.read_text() == added

    assert not update_release_history(
        Version("6.8.0"), date(2026, 9, 22), path=release_notes
    )
    assert release_notes.read_text() == added


def test_get_release_history_doi(release_notes):
    assert (
        get_release_history_doi(Version("6.7.0"), path=release_notes)
        == "https://doi.org/10.5066/P1IJAXDZ"
    )
    assert (
        get_release_history_doi(Version("6.8.0"), path=release_notes)
        == "https://doi.org/10.5066/P1PGE9XW"
    )
    assert get_release_history_doi(Version("6.8.1"), path=release_notes) is None
    assert get_release_history_doi(Version("6.9.0.dev0"), path=release_notes) is None

    # a row added for a patch release is found, with the reused DOI
    update_release_history(Version("6.8.1"), date(2026, 9, 21), path=release_notes)
    assert (
        get_release_history_doi(Version("6.8.1"), path=release_notes)
        == "https://doi.org/10.5066/P1PGE9XW"
    )


def test_update_release_history_no_table(tmp_path):
    path = tmp_path / "ReleaseNotes.tex"
    path.write_text("nothing to see here\n")
    with pytest.raises(ValueError, match="No release history table rows"):
        update_release_history(Version("6.8.1"), date(2026, 9, 21), path=path)


def test_update_release_history_repo_notes():
    # the table in the repository's release notes must be parseable
    assert list(_release_history_row.finditer(release_notes_path.read_text()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
Update version information stored in version.txt in the project root,
as well as several other files in the repository:

  ../version.txt
  ../meson.build
  ../utils/mf5to6/meson.build
  ../doc/version.tex
  ../README.md
  ../DISCLAIMER.md
  ../code.json
  ../src/Utilities/version.f90.in
  ../src/Utilities/version.f90

These include a combination of version strings, build timestamps, disclaimer
text, text indicating whether the release is provisional or approved, source
code setting the variable IDEVELOPMODE to either 0 or 1, and other data.

Provide a `--version` string following semantic versioning conventions.
If none of --version, --release or --post-release is provided, the version
number will not be changed, just timestamps.

Use `--release` (`-r`) to use the current development version with any
development segment (e.g. '.dev0') removed. Use `--post-release` (`-p`) to
use the next development version, with the minor version incremented and a
'.dev0' suffix re-added; this is the version the post-release reset sets on
the develop branch.

Use `--get` (`-g`) or `--dry-run` to print the resolved version without
making changes. The resolved version is printed on the last line in every
mode, so `--post-release` alone both updates the files and reports the new
version.

Use `--releasemode` to control whether IDEVELOPMODE is set to 0 instead
of 1, and to alter mf6's output and disclaimer text reflecting approval.

Use `--citation` (`-c`) to render the current software citation. Pass the
release DOI link via `--doi` (`-d`), e.g.
`--doi https://doi.org/10.5066/P1PGE9XW`; if omitted, the DOI in the release
history row for the version in ReleaseNotes.tex is used, or if there is no
row, the umbrella MODFLOW software DOI.

When updating to a release version (not a pre-release or development version),
a row for the version is added to the release history table in ReleaseNotes.tex
unless it already has one. The row is dated today, or `--date` (`YYYY-MM-DD`)
if provided. A patch release uses the DOI of the previous release with the same
major and minor version. A new minor or major release needs `--doi`, unless its
row was added to ReleaseNotes.tex beforehand.
            """
        ),
    )
    parser.add_argument(
        "-c",
        "--citation",
        required=False,
        action="store_true",
        help="Show the citation, don't update anything. Defaults to false.",
    )
    parser.add_argument(
        "-g",
        "--get",
        required=False,
        action="store_true",
        help="Show the version, don't update anything. Defaults to false",
    )
    parser.add_argument(
        "-d",
        "--doi",
        required=False,
        default=None,
        help="DOI link (e.g. https://doi.org/10.5066/P1PGE9XW) of the release. "
        "Substituted into the software citation rendered by --citation, "
        "defaulting to the DOI in the release history row for the version, or "
        f"if there is none the umbrella MODFLOW software DOI ({_default_doi}). "
        "Also used in the release history row added when updating to a release "
        "version, defaulting to the DOI of the previous release with the same "
        "major and minor version. Required for a new minor or major release "
        "whose release history row doesn't already exist.",
    )
    parser.add_argument(
        "--date",
        required=False,
        type=date.fromisoformat,
        default=None,
        help="Release date (YYYY-MM-DD) for the release history row added when "
        "updating to a release version. Defaults to today.",
    )
    parser.add_argument(
        "-a",
        "--releasemode",
        required=False,
        action="store_true",
        help="Enable release mode for a full, standard release. Modifies "
        "disclaimer language reflecting approval. Sets IDEVELOPMODE = 0. "
        "Defaults to false for preliminary development distributions.",
    )
    parser.add_argument(
        "-v",
        "--version",
        required=False,
        help="Specify the release version. Value must follow PEP 440.",
    )
    parser.add_argument(
        "-r",
        "--release",
        required=False,
        action="store_true",
        help="Use the current development version with any development segment "
        "(e.g. '.dev0') removed. Defaults to false.",
    )
    parser.add_argument(
        "-p",
        "--post-release",
        required=False,
        action="store_true",
        help="Use the next development version, with the minor version "
        "incremented and a '.dev0' suffix re-added. Defaults to false.",
    )
    parser.add_argument(
        "--dry-run",
        required=False,
        action="store_true",
        help="Print the resolved version and exit without updating anything. "
        "Defaults to false.",
    )

    args = parser.parse_args()
    citation = args.citation
    developmode = not args.releasemode

    if args.post_release:
        version = post_release_version()
    elif args.release:
        version = release_version()
    elif args.version:
        version = Version(args.version)
    else:
        version = _current_version

    if citation:
        doi = args.doi or get_release_history_doi(version)
        if doi is None:
            print(
                f"No release history row for {version}, using the umbrella DOI",
                file=sys.stderr,
            )
        print(
            get_software_citation(
                timestamp=datetime.now(),
                version=version,
                doi=doi or _default_doi,
                developmode=developmode,
            )
        )
    elif args.get or args.dry_run:
        print(version)
    else:
        mode = "develop" if developmode else "release"
        print(f"Updating to version {version} in {mode} mode", file=sys.stderr)
        update_version(
            version=version,
            timestamp=datetime.now(),
            developmode=developmode,
            release_date=args.date,
            doi=args.doi,
        )
        print(version)

#!/usr/bin/python

"""
Update files in this modflow6 repository according to release information.

This script is used to update several files in the modflow6 repository, including:

  ../version.txt
  ../meson.build
  ../doc/version.tex
  ../README.md
  ../DISCLAIMER.md
  ../code.json
  ../src/Utilities/version.f90

Information in these files include version number (major.minor.patch[label]), build
timestamp, whether or not the release is preliminary/provisional or official/approved,
whether the source code should be compiled in develop mode (IDEVELOPMODE = 1) or for
release, and other metadata.

The version number is read from ../version.txt, which contains major, minor, and patch
version numbers, and an optional label. Version numbers are substituted into source
code, latex files, markdown files, etc. The version number can be provided explicitly
using --version, short -v.

If --releasemode is provided, IDEVELOPMODE is set to 0 in src/Utilities/version.f90.
Otherwise, IDEVELOPMODE is set to 1.

if --releasemode is provided, the disclaimer in src/Utilities/version.f90 and the
README/DISCLAIMER markdown files is modified to reflect review and approval.
Otherwise the language reflects preliminary/provisional
status, and version strings contain "(preliminary)".
"""

import argparse
import json
import os
import textwrap
from collections import OrderedDict
from datetime import datetime
from pathlib import Path
from typing import Optional

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
    project_root_path / "doc" / "version.tex",
    project_root_path / "doc" / "version.py",
    project_root_path / "README.md",
    project_root_path / "DISCLAIMER.md",
    project_root_path / "CITATION.cff",
    project_root_path / "code.json",
    project_root_path / "src" / "Utilities" / "version.f90",
]


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


def get_software_citation(
    timestamp: datetime, version: Version, developmode: bool = False
) -> str:
    # get data Software/Code citation for FloPy
    citation = yaml.safe_load((project_root_path / "CITATION.cff").read_text())

    sb = ""
    if developmode:
        sb = " (preliminary)"
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
        f"MODFLOW 6 Modular Hydrologic Model version {version}{sb}: "
        f"U.S. Geological Survey Software Release, {timestamp:%-d %B %Y}, "
        "https://doi.org/10.5066/P9FL1JCC"
    )

    return line


def log_update(path, version: Version):
    print(f"Updated {path} with version {version}")


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
    path = project_root_path / "meson.build"
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
    version: Optional[Version],
    timestamp: datetime,
    developmode: bool = False,
):
    path = project_root_path / "src" / "Utilities" / "version.f90"
    lines = open(path, "r").read().splitlines()
    with open(path, "w") as f:
        skip = False
        version_spl = str(version).rpartition("-")
        if version_spl[1]:
            version_num = version_spl[0]
            version_label = version_spl[2]
        else:
            version_num = str(version_spl[2])
            version_label = ""

        for line in lines:
            # skip all of the disclaimer text
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
            elif ":: VERSIONTAG =" in line:
                fmat_tstmp = timestamp.strftime("%m/%d/%Y")
                label_clause = version_label if version_label else ""
                label_clause += " (preliminary)" if developmode else ""
                line = (
                    line.rpartition("::")[0]
                    + f":: VERSIONTAG = '{label_clause} {fmat_tstmp}'"
                )
            elif ":: FMTDISCLAIMER =" in line:
                line = get_disclaimer(developmode=developmode, formatted=True)
                skip = True
            f.write(f"{line}\n")
    log_update(path, version)


def update_readme_and_disclaimer(version: Version, developmode: bool = False):
    disclaimer = get_disclaimer(developmode, formatted=False)
    readme_path = str(project_root_path / "README.md")
    readme_lines = open(readme_path, "r").read().splitlines()
    with open(readme_path, "w") as f:
        for line in readme_lines:
            if "## Version " in line:
                version_line = f"### Version {version}"
                if developmode:
                    version_line += " (preliminary)"
                f.write(f"{version_line}\n")
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


def update_version(
    version: Version = None,
    timestamp: datetime = datetime.now(),
    developmode: bool = False,
):
    """
    Update version information stored in version.txt in the project root,
    as well as several other files in the repository. Version updates are
    performed by explicitly providing a version argument to this function
    and a lock is held on the version file to make sure that the state of
    the multiple files containing version information stays synchronized.
    If no version argument is provided, the version number isn't changed.
    """

    lock_path = Path(version_file_path.name + ".lock")
    try:
        lock = FileLock(lock_path)
        previous = Version(version_file_path.read_text().strip())
        version = version if version else previous

        with lock:
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
        update_version(timestamp=timestamp, version=version, developmode=full)
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
        version_f90_path = project_root_path / "src" / "Utilities" / "version.f90"
        lines = version_f90_path.read_text().splitlines()
        assert any(f"IDEVELOPMODE = {0 if full else 1}" in line for line in lines)

        # check disclaimer has appropriate language
        disclaimer_path = project_root_path / "DISCLAIMER.md"
        lines = disclaimer_path.read_text().splitlines()
        assert any(("approved for release") in line for line in lines) == full
        assert any(("preliminary or provisional") in line for line in lines) != full

        # check readme has appropriate language
        readme_path = project_root_path / "README.md"
        lines = readme_path.read_text().splitlines()
        assert any(("(preliminary)") in line for line in lines) != full
    finally:
        for p in touched_file_paths:
            os.system(f"git restore {p}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
Update version information stored in version.txt in the project root,
as well as several other files in the repository:

  ../version.txt
  ../meson.build
  ../doc/version.tex
  ../README.md
  ../DISCLAIMER.md
  ../code.json
  ../src/Utilities/version.f90

These include a combination of version strings, build timestamps, disclaimer
text, text indicating whether the release is provisional or approved, source
code setting the variable IDEVELOPMODE to either 0 or 1, and other data.

Provide a `--version` string following semantic versioning conventions.
If --version is not provided, the version number will not be changed,
just timestamps.

Use `--get` (`-g`) to show the current version without making changes.
The version number is read from version.txt in the project root.

Use `--releasemode` to control whether IDEVELOPMODE is set to 0 instead
of 1, and to alter mf6's output and disclaimer text reflecting approval.

Use `--citation` (`-c`) to render the current software citation.
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

    args = parser.parse_args()
    get = args.get
    citation = args.citation
    developmode = not args.releasemode
    version = Version(args.version) if args.version else _current_version

    if get:
        print(Version((project_root_path / "version.txt").read_text().strip()))
    elif citation:
        print(
            get_software_citation(
                timestamp=datetime.now(), version=version, developmode=developmode
            )
        )
    else:
        mode = "develop" if developmode else "release"
        print(f"Updating to version {version} in {mode} mode")
        update_version(
            version=version, timestamp=datetime.now(), developmode=developmode
        )

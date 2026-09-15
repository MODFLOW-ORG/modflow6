"""Post-release housekeeping for the release notes.

Archives the just-released version's notes to previous/v<version>.tex in the
compact archive format, registers them in appendixA.tex, and clears the
archived items from develop.toml so the next development cycle starts empty.

The version is read from version.txt, so run this before bumping the version.
For a patch release (--patch) only the fix items are archived and cleared; the
rest carry forward to the next minor release.
"""

import argparse
import re
import sys
from warnings import warn

from mk_releasenotes import latest_release, notes_dir, render, version

develop_toml_path = notes_dir / "develop.toml"
appendix_path = notes_dir / "appendixA.tex"
previous_dir = notes_dir / "previous"


def register_archive(version: str):
    """Add an \\input line for previous/v<version>.tex to appendixA.tex, just
    above the newest existing entry. No-op if the line is already present."""
    input_line = f"\\input{{./previous/v{version}.tex}}"
    lines = appendix_path.read_text().splitlines()
    if any(input_line in line for line in lines):
        return
    for i, line in enumerate(lines):
        if "\\input{./previous/" in line:
            indent = line[: len(line) - len(line.lstrip())]
            lines.insert(i, f"{indent}{input_line}")
            break
    else:
        lines.append(f"        {input_line}")
    appendix_path.write_text("\n".join(lines) + "\n")
    print(f"Registered {input_line} in {appendix_path}", file=sys.stderr)


def clear_develop_toml(*, patch: bool = False):
    """Remove release note items from develop.toml for the next cycle.

    The [sections] and [subsections] tables are kept. For a patch release,
    non-fix items are kept too (they carry forward to the next minor release);
    otherwise every item is removed.
    """
    lines = develop_toml_path.read_text().splitlines()
    try:
        first = next(i for i, ln in enumerate(lines) if ln.strip() == "[[items]]")
    except StopIteration:
        return  # already empty

    header = lines[:first]
    while header and not header[-1].strip():
        header.pop()

    blocks: list[list[str]] = []
    for ln in lines[first:]:
        if ln.strip() == "[[items]]":
            blocks.append([ln])
        else:
            blocks[-1].append(ln)

    kept = []
    if patch:
        for block in blocks:
            section = next(
                (
                    m.group(1)
                    for ln in block
                    if (m := re.match(r'\s*section\s*=\s*"([^"]*)"', ln))
                ),
                "",
            )
            if section != "fixes":
                while block and not block[-1].strip():
                    block.pop()
                kept.append(block)

    text = "\n".join(header) + "\n"
    for block in kept:
        text += "\n" + "\n".join(block) + "\n"
    develop_toml_path.write_text(text)
    print(
        f"Cleared {develop_toml_path}: removed {len(blocks) - len(kept)} item(s), "
        f"kept {len(kept)}",
        file=sys.stderr,
    )


def reset_release_notes(*, patch: bool = False):
    """Archive the just-released version's notes and clear develop.toml."""
    header_version, header_date = latest_release()
    if header_version != version:
        warn(
            f"version.txt is {version} but the last row of the ReleaseNotes.tex "
            f"releases table is {header_version}; archiving as {version}. Add the "
            "release row before releasing (see the release procedure)."
        )

    if render(
        develop_toml_path,
        previous_dir / f"v{version}.tex",
        patch=patch,
        archive=True,
        version=version,
        date=header_date,
    ):
        print(
            f"Archived release notes to {previous_dir / f'v{version}.tex'}",
            file=sys.stderr,
        )
        register_archive(version)
    else:
        warn(f"No {'fix ' if patch else ''}items to archive for v{version}")

    clear_develop_toml(patch=patch)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--patch",
        default=False,
        action="store_true",
        help="Archive and clear only the fix items (patch release).",
    )
    args = parser.parse_args()
    reset_release_notes(patch=args.patch)

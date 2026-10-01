"""Post-release housekeeping for the release notes.

Archives the just-released version's notes to previous/v<version>.tex in the
compact archive format, registers them in appendixA.tex, and deletes the
archived item files from items/ so the next development cycle starts empty.

The version is read from version.txt, so run this before bumping the version.
For a patch release (--patch) only the fix and example items are archived and
cleared; the rest carry forward to the next minor release.
"""

import argparse
import sys
from warnings import warn

from mk_releasenotes import (
    items_dir,
    latest_release,
    load_items,
    notes_dir,
    patch_sections,
    render,
    tomllib,
    version,
)

schema_path = notes_dir / "schema.toml"
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


def clear_items(*, patch: bool = False):
    """Delete release note item files from items/ for the next cycle.

    For a patch release, items outside the patch sections (fixes and examples)
    are kept (they carry forward to the next minor release); otherwise every
    item is deleted. The README.md is kept so the directory stays tracked.
    """
    with open(schema_path, "rb") as schema_file:
        schema = tomllib.load(schema_file)
    items = load_items(items_dir, schema["sections"], schema["subsections"])
    removed = [
        path for path, item in items if not patch or item["section"] in patch_sections
    ]
    for path in removed:
        path.unlink()
    print(
        f"Cleared {items_dir}: removed {len(removed)} item(s), "
        f"kept {len(items) - len(removed)}",
        file=sys.stderr,
    )


def reset_release_notes(*, patch: bool = False):
    """Archive the just-released version's notes and clear the items."""
    header_version, header_date = latest_release()
    if header_version != version:
        warn(
            f"version.txt is {version} but the last row of the ReleaseNotes.tex "
            f"releases table is {header_version}; archiving as {version}. Add the "
            "release row before releasing (see the release procedure)."
        )

    if render(
        schema_path,
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
        warn(f"No {'patch ' if patch else ''}items to archive for v{version}")

    clear_items(patch=patch)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--patch",
        default=False,
        action="store_true",
        help="Archive and clear only the fix and example items (patch release).",
    )
    args = parser.parse_args()
    reset_release_notes(patch=args.patch)

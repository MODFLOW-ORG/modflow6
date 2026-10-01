"""Convert the release note items to a LaTeX file for the PDF build.

Each release note item is a TOML file in items/ with a section, subsection
and description. Valid sections and subsections are defined in schema.toml,
in the order they are rendered.

Two formats (see --archive). The --archive format is more compact, for the
archive section of the release notes document, and has a leading version string
header (read from the last row of the releases table in ReleaseNotes.tex). The
default format, for the release notes document section, omits that header and
uses more widely spaced section headers.

See reset_releasenotes.py for the post-release archive-and-clear step.
"""

import argparse
import datetime
import sys
from pathlib import Path
from warnings import warn

import release_history

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

notes_dir = Path(__file__).parent
items_dir = notes_dir / "items"
version_file = Path(__file__).parents[2] / "version.txt"
version = version_file.read_text().strip()
date = datetime.date.today().strftime("%b %d, %Y")

# sections included in the release notes for a patch release. Bug fixes always
# ship in a patch. New examples are included too, since examples are versioned
# separately from the program and have shipped in patch releases before.
patch_sections = ("fixes", "examples")


def latest_release():
    """Version and date of the most recent release, read from the last row
    of the releases table in ReleaseNotes.tex."""
    try:
        return release_history.latest_release()
    except ValueError as e:
        raise ValueError(f"{e}; pass --version and --date explicitly") from e


def load_items(
    items_dir: Path, sections: dict, subsections: dict
) -> list[tuple[Path, dict]]:
    """Load and validate the release note item files in a directory.

    Returns (path, item) pairs sorted by file name. Raises ValueError listing
    every invalid item if any is missing a required key or has a section or
    subsection not defined in the schema.
    """
    items = []
    errors = []
    for path in sorted(items_dir.glob("*.toml")):
        with open(path, "rb") as item_file:
            item = tomllib.load(item_file)
        for key in ("section", "subsection", "description"):
            if key not in item:
                errors.append(f"{path.name}: missing required key '{key}'")
        if "section" in item and item["section"] not in sections:
            errors.append(
                f"{path.name}: invalid section '{item['section']}'"
                f", expected one of: {list(sections)}"
            )
        if item.get("subsection") and item["subsection"] not in subsections:
            errors.append(
                f"{path.name}: invalid subsection '{item['subsection']}'"
                f", expected one of: {list(subsections)}"
            )
        items.append((path, item))
    if errors:
        raise ValueError("Invalid release note items:\n" + "\n".join(errors))
    return items


def render(
    toml_path: Path,
    tex_path: Path,
    *,
    template_name: str = "develop.tex.jinja",
    patch: bool = False,
    archive: bool = False,
    version: str = version,
    date: str = date,
) -> bool:
    """Render the release note items to a LaTeX file, using the sections and
    subsections in the schema TOML file at toml_path. Items are read from the
    items/ directory next to the schema file.

    Returns True if notes were rendered, False if there was nothing to render
    (no schema file, or no items after any --patch filtering). In the latter case
    an empty LaTeX file is still written so downstream document builds succeed.
    """
    if not toml_path.is_file():
        warn(f"Release notes schema file not found: {toml_path}")
        return False

    tex_path.unlink(missing_ok=True)

    from jinja2 import Environment, FileSystemLoader

    with open(toml_path, "rb") as toml_file:
        content = tomllib.load(toml_file)
    sections = content.get("sections", {})
    subsections = content.get("subsections", {})
    items = [
        item
        for _, item in load_items(toml_path.parent / "items", sections, subsections)
    ]
    # if patch, only include fixes and examples
    if patch:
        items = [item for item in items if item["section"] in patch_sections]
        sections = {k: v for k, v in sections.items() if k in patch_sections}
        used = {item.get("subsection") for item in items}
        subsections = {k: v for k, v in subsections.items() if k in used}
    # make sure each item has a subsection entry even if empty
    for item in items:
        if not item.get("subsection"):
            item["subsection"] = ""
    # items without a subsection come first in their section, with no header
    subsections = {"": "", **subsections}
    if not any(items):
        warn("No release notes found, aborting")
        # still leave an empty file behind
        tex_path.write_text("")
        return False

    loader = FileSystemLoader(notes_dir)
    env = Environment(
        loader=loader,
        trim_blocks=True,
        lstrip_blocks=True,
        line_statement_prefix="_",
        keep_trailing_newline=False,
        # since latex uses curly brackets,
        # replace block/var start/end tags
        block_start_string="([",
        block_end_string="])",
        variable_start_string="((",
        variable_end_string="))",
    )
    template = env.get_template(template_name)
    rendered = template.render(
        sections=sections,
        subsections=subsections,
        items=items,
        version=version,
        date=date,
        archive=archive,
    )
    tex_path.write_text(rendered.rstrip() + "\n")
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--toml", default="schema.toml")
    parser.add_argument("--tex", default="develop.tex")
    parser.add_argument("--patch", default=False, action="store_true")
    parser.add_argument(
        "--archive",
        default=False,
        action="store_true",
        help=(
            "Render this version's release notes in a more compact format for the "
            "archive section of the release notes document. Version/date are read "
            "from the last row of the releases table in ReleaseNotes.tex unless a "
            "--version and/or --date are given. The default rendering format, for "
            "the release notes document section, omits the leading version string "
            "header and uses more widely spaced section headers."
        ),
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Override the version string in the --archive header.",
    )
    parser.add_argument(
        "--date",
        default=None,
        help="Override the date in the --archive header.",
    )
    args = parser.parse_args()

    toml_path = Path(args.toml).expanduser().absolute()
    tex_path = Path(args.tex).expanduser().absolute()

    render_version = version
    render_date = date
    if args.archive:
        release_version, release_date = latest_release()
        render_version = args.version or release_version
        render_date = args.date or release_date

    try:
        render(
            toml_path,
            tex_path,
            template_name=f"{tex_path.name}.jinja",
            patch=args.patch,
            archive=args.archive,
            version=render_version,
            date=render_date,
        )
    except ValueError as e:
        sys.exit(str(e))

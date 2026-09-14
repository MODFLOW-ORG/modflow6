"""Convert the release notes TOML file to a LaTeX file for the PDF build.

Two formats (see --archive). The --archive format is more compact, for the
archive section of the release notes document, and has a leading version string
header (read from the last row of the releases table in ReleaseNotes.tex). The
default format, for the release notes document section, omits that header and
uses more widely spaced section headers.
"""

import argparse
import datetime
import re
import sys
from pathlib import Path
from warnings import warn

notes_dir = Path(__file__).parent
version_file = Path(__file__).parents[2] / "version.txt"
version = version_file.read_text().strip()
date = datetime.date.today().strftime("%b %d, %Y")


def latest_release():
    """Version and date of the most recent release, read from the last row
    of the releases table in ReleaseNotes.tex."""
    rows = re.findall(
        r"^\s*(\d+\.\d+\.\d+)\s*&\s*([^&]+?)\s*&\s*\\url",
        (notes_dir / "ReleaseNotes.tex").read_text(),
        re.MULTILINE,
    )
    if not rows:
        raise ValueError(
            "No rows found in the releases table in ReleaseNotes.tex; "
            "pass --version and --date explicitly"
        )
    return rows[-1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
    patch = args.patch
    archive = args.archive
    if archive:
        release_version, release_date = latest_release()
        version = args.version or release_version
        date = args.date or release_date
    if not toml_path.is_file():
        warn(f"Release notes TOML file not found: {toml_path}")
        sys.exit(0)

    tex_path.unlink(missing_ok=True)

    import tomli
    from jinja2 import Environment, FileSystemLoader

    with open(toml_path, "rb") as f:
        content = tomli.load(f)
    sections = content.get("sections", {})
    subsections = content.get("subsections", {})

    items_dir = toml_path.parent / "items"
    items = []
    errors = []
    for fragment in sorted(items_dir.glob("*.toml")):
        with open(fragment, "rb") as f:
            item = tomli.load(f)
        for key in ("section", "subsection", "description"):
            if key not in item:
                errors.append(f"{fragment.name}: missing required key '{key}'")
        if "section" in item and item["section"] not in sections:
            errors.append(
                f"{fragment.name}: invalid section '{item['section']}'"
                f", expected one of: {list(sections)}"
            )
        if (
            "subsection" in item
            and item["subsection"]
            and item["subsection"] not in subsections
        ):
            errors.append(
                f"{fragment.name}: invalid subsection '{item['subsection']}'"
                f", expected one of: {list(subsections)}"
            )
        items.append(item)
    if errors:
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)

    # if patch, only include fixes
    if patch:
        items = [item for item in items if item["section"] == "fixes"]
        sections = {k: v for k, v in sections.items() if k == "fixes"}
        used_subsections = {item["subsection"] for item in items}
        subsections = {k: v for k, v in subsections.items() if k in used_subsections}
    # make sure each item has a subsection entry even if empty
    for item in items:
        if not item.get("subsection"):
            item["subsection"] = ""
    if not any(items):
        warn("No release notes found, aborting")
        # still leave an empty file behind so downstream document builds succeed
        tex_path.write_text("")
        sys.exit(0)

    loader = FileSystemLoader(Path(__file__).parent)
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
    template = env.get_template(f"{tex_path.name}.jinja")
    with open(tex_path, "w") as tex_file:
        rendered = template.render(
            sections=sections,
            subsections=subsections,
            items=items,
            version=version,
            date=date,
            archive=archive,
        )
        tex_file.write(rendered.rstrip() + "\n")

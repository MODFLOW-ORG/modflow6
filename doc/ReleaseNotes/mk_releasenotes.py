# This script converts the release notes TOML file
# to a latex file, from which is later built a PDF.
import argparse
import datetime
import sys
from pathlib import Path
from warnings import warn

version_file = Path(__file__).parents[2] / "version.txt"
version = version_file.read_text().strip()
date = datetime.date.today().strftime("%b %d, %Y")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--toml", default="schema.toml")
    parser.add_argument("--tex", default="develop.tex")
    parser.add_argument("--patch", default=False, action="store_true")
    args = parser.parse_args()
    toml_path = Path(args.toml).expanduser().absolute()
    tex_path = Path(args.tex).expanduser().absolute()
    patch = args.patch
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

    changes_dir = toml_path.parent / "changes"
    items = []
    errors = []
    for fragment in sorted(changes_dir.glob("*.toml")):
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
        tex_file.write(
            template.render(
                sections=sections,
                subsections=subsections,
                items=items,
                version=version,
                date=date,
            )
        )

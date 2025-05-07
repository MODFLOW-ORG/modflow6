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
    parser.add_argument("--toml", default="develop.toml")
    parser.add_argument("--tex", default="develop.tex")
    args = parser.parse_args()
    toml_path = Path(args.toml).expanduser().absolute()
    tex_path = Path(args.tex).expanduser().absolute()
    if not toml_path.is_file():
        warn(f"Release notes TOML file not found: {toml_path}")
        sys.exit(0)

    tex_path.unlink(missing_ok=True)

    import tomli
    from jinja2 import Environment, FileSystemLoader

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
        with open(toml_path, "rb") as toml_file:
            content = tomli.load(toml_file)
            sections = content.get("sections", [])
            subsections = content.get("subsections", [])
            items = content.get("items", [])
            # make sure each item has a subsection entry even if empty
            for item in items:
                if not item.get("subsection"):
                    item["subsection"] = ""
            if not any(items):
                warn("No release notes found, aborting")
                sys.exit(0)
            tex_file.write(
                template.render(
                    sections=sections,
                    subsections=subsections,
                    items=items,
                    version=version,
                    date=date,
                )
            )

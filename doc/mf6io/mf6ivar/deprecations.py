from pathlib import Path
from typing import Optional

from packaging.version import Version

PROJ_ROOT_PATH = Path(__file__).parents[3]
MF6IVAR_PATH = PROJ_ROOT_PATH / "doc" / "mf6io" / "mf6ivar"


def get_deprecations(
    dfndir,
) -> list[tuple[Path, str, Version, Optional[Version]]]:
    dfns = Path(dfndir).rglob("*.dfn")
    deps = {}
    for dfn in dfns:
        with open(dfn, "r") as f:
            name = None
            for line in f:
                if line.startswith("#"):
                    continue
                keys = ["deprecated", "removed"]
                ikeys = {k: i for i, k in enumerate(keys)}
                for key in keys:
                    if line.startswith("name"):
                        name = line.split()[1]
                    if line.startswith(key):
                        val = deps.get((dfn, key), [None, None])
                        key, ver = line.split()
                        ik = ikeys[key]
                        val[ik] = val[ik] if val[ik] else Version(ver)
                        deps[(dfn, name)] = val

    return [(file, key, dep, rem) for (file, key), (dep, rem) in deps.items()]


def create_deprecations_file(dfndir, mddir, verbose):
    deprecations = get_deprecations(dfndir)
    deps_path = (mddir / "deprecations.md").absolute()
    if verbose:
        print(f"Found {len(deprecations)} deprecations, writing {deps_path}")
    with open(deps_path, "w") as f:
        s = "#### Deprecations\n\n"
        s += (
            "The following table lists deprecated options and the versions "
            "in which they were deprecated and (optionally) removed.\n\n"
        )
        if any(deprecations):
            s += "| Model-Package | Option | Deprecated | Removed |\n"
            s += "|:--------------|:-------|:-----------|:--------|\n"
            for file, option, deprecated, removed in deprecations:
                deprecated = deprecated if deprecated else removed if removed else ""
                removed = removed if removed else ""
                s += f"| {file.stem} | {option} | {deprecated} | {removed} |\n"
            if len(s) > 0:
                s += "\n"
        f.write(s)


if __name__ == "__main__":
    dfndir = MF6IVAR_PATH / "dfn"
    mddir = MF6IVAR_PATH / "md"
    mddir.mkdir(exist_ok=True)
    create_deprecations_file(dfndir, mddir, verbose=True)

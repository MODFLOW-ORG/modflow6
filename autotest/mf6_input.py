"""
Reading MODFLOW 6 input files and DFN definitions as text, for the re-run
comparisons that rewrite a test simulation's input.

Input files are read line by line: comment lines (starting with `#` or
`!`) are skipped, and a block is found by its BEGIN and END lines. Inline
comments and quoted names containing spaces aren't handled.
"""

from pathlib import Path

DFN_PATH = Path(__file__).parents[1] / "doc" / "mf6io" / "mf6ivar" / "dfn"


def read_lines(path: Path) -> list[str]:
    """The lines of an input file."""
    return Path(path).read_text().splitlines()


def tokens(line: str) -> list[str]:
    """The tokens of an input file line, none for a comment line."""
    words = line.split()
    return [] if not words or words[0].startswith(("#", "!")) else words


def block_span(lines: list[str], name: str) -> tuple[int, int] | None:
    """Line indices of the BEGIN and END lines of the first block `name`,
    or None if there is no such block."""
    name = name.upper()
    begin = None
    for i, line in enumerate(lines):
        upper = [t.upper() for t in tokens(line)[:2]]
        if upper == ["BEGIN", name] and begin is None:
            begin = i
        elif upper == ["END", name] and begin is not None:
            return begin, i
    return None


def block_rows(lines: list[str], name: str, source: Path | str = "") -> list[list[str]]:
    """The tokens of each line of the first block `name`, without empty
    and comment lines. `source` names the file in the error raised if
    there is no such block."""
    span = block_span(lines, name)
    if span is None:
        raise ValueError(f"No {name.upper()} block in {source or 'input'}")
    rows = (tokens(line) for line in lines[span[0] + 1 : span[1]])
    return [row for row in rows if row]


def models(ws: Path) -> list[tuple[str, str, str]]:
    """(model type, name file, model name) of each model in a simulation
    workspace, the model type without its version, e.g. "gwf"."""
    path = Path(ws) / "mfsim.nam"
    return [
        (r[0].lower().removesuffix("6"), r[1], r[2])
        for r in block_rows(read_lines(path), "MODELS", path)
    ]


def packages(ws: Path, namefile: str, mtype: str) -> list[tuple[str, str, str]]:
    """(package type, file, package name) of each package of a model. A
    package the name file gives no name is named as MODFLOW 6 names it: by
    its type, or for a type a model can have several of, `<type>-<n>`, n
    counting the packages of that type."""
    ws = Path(ws)
    path = ws / namefile
    result, count = [], {}
    for row in block_rows(read_lines(path), "PACKAGES", path):
        ptype = row[0].lower().removesuffix("6")
        count[ptype] = count.get(ptype, 0) + 1
        if len(row) > 2:
            pname = row[2].lower()
        elif multi_package(component(mtype, ptype, read_lines(ws / row[1]))):
            pname = f"{ptype}-{count[ptype]}"
        else:
            pname = ptype
        result.append((ptype, row[1], pname))
    return result


def component(mtype: str, ptype: str, lines: list[str]) -> str | None:
    """The DFN component of a package file, e.g. "gwf-rcha" for a GWF RCH
    package reading its input as arrays, or None if there is none."""
    words = {t.upper() for line in lines for t in tokens(line)}
    variant = (
        "a" if "READASARRAYS" in words else "g" if "READARRAYGRID" in words else ""
    )
    name = f"{mtype}-{ptype}{variant}"
    return name if (DFN_PATH / f"{name}.dfn").is_file() else None


def multi_package(component: str | None) -> bool:
    """Whether a model can have several packages of a DFN component, as
    the DFN's `flopy multi-package` header comment says."""
    if component is None:
        return False
    head = (DFN_PATH / f"{component}.dfn").read_text().splitlines()[:5]
    return any("multi-package" in line for line in head)


def dfn_names(component: str, block: str) -> set[str]:
    """Names of the variables in a block of a DFN component."""
    names, current = set(), None
    for line in (DFN_PATH / f"{component}.dfn").read_text().splitlines():
        words = line.split()
        if words[:1] == ["block"] and len(words) > 1:
            current = words[1]
        elif words[:1] == ["name"] and len(words) > 1 and current == block:
            names.add(words[1])
    return names


def netcdf_arrays(component: str) -> set[str]:
    """Names of the arrays of a DFN component MODFLOW 6 can read from
    NetCDF."""
    names, name = set(), None
    for line in (DFN_PATH / f"{component}.dfn").read_text().splitlines():
        words = line.split()
        if words[:1] == ["name"] and len(words) > 1:
            name = words[1]
        elif words == ["netcdf", "true"] and name:
            names.add(name)
    return names


def aux_names(lines: list[str]) -> list[str]:
    """A package's auxiliary variable names, in order."""
    span = block_span(lines, "OPTIONS")
    for line in lines[span[0] + 1 : span[1]] if span else []:
        words = tokens(line)
        if words and words[0].upper() in ("AUXILIARY", "AUX"):
            return [t.lower() for t in words[1:]]
    return []


def fmi_input_files(ws: Path) -> set[str]:
    """Names of the files the FMI packages of a simulation's models read,
    input even when named like output (e.g. a GWF budget file)."""
    ws = Path(ws)
    if not (ws / "mfsim.nam").is_file():
        return set()
    files = set()
    for mtype, namefile, _ in models(ws):
        for ptype, pfile, _ in packages(ws, namefile, mtype):
            lines = read_lines(ws / pfile) if ptype == "fmi" else []
            # an FMI package of a model coupled to its flow model by an
            # exchange reads no files
            if block_span(lines, "PACKAGEDATA") is None:
                continue
            for row in block_rows(lines, "PACKAGEDATA"):
                upper = [t.upper() for t in row]
                if "FILEIN" in upper[:-1]:
                    files.add(Path(row[upper.index("FILEIN") + 1].strip("\"'")).name)
    return files

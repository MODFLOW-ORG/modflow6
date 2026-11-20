"""
Package input read from the NetCDF file MODFLOW 6 itself exports.

For `Comparison.MF6_EXPORT_STRUCTURED` and `Comparison.MF6_EXPORT_LAYERED`:
a copy of a test simulation is run in validate mode with every package that
can export its arrays doing so (`enable_export`), then its package files are
changed to read those arrays from the exported file with the NETCDF keyword
(`use_netcdf_input`), and it is run again. This tests NetCDF input without
flopy4, which the FP4 comparisons need and the MODFLOW 6 test runs don't
install.

Arrays are replaced where they are given, in the GRIDDATA block and in each
PERIOD block, so the PERIOD blocks keep their periods and the order of their
arrays. TIMEARRAYSERIES entries are kept.
"""

from pathlib import Path

DFN_PATH = Path(__file__).parents[1] / "doc" / "mf6io" / "mf6ivar" / "dfn"
OUTPUT_OPTION = {"structured": "NETCDF_STRUCTURED", "layered": "NETCDF_MESH2D"}


def _lines(path: Path) -> list[str]:
    return path.read_text().splitlines()


def _tokens(line: str) -> list[str]:
    tokens = line.split()
    return [] if not tokens or tokens[0].startswith(("#", "!")) else tokens


def _block(lines: list[str], name: str) -> tuple[int, int] | None:
    """Line indices of the BEGIN and END lines of the first block `name`."""
    begin = None
    for i, line in enumerate(lines):
        upper = [t.upper() for t in _tokens(line)[:2]]
        if upper == ["BEGIN", name] and begin is None:
            begin = i
        elif upper == ["END", name] and begin is not None:
            return begin, i
    return None


def _dfn_names(component: str, block: str) -> set[str]:
    """Names of the variables in a block of a DFN component."""
    names, name, current = set(), None, None
    for row in (DFN_PATH / f"{component}.dfn").read_text().splitlines():
        tokens = row.split()
        if tokens[:1] == ["block"] and len(tokens) > 1:
            current = tokens[1]
        elif tokens[:1] == ["name"] and len(tokens) > 1:
            name = tokens[1]
            if current == block:
                names.add(name)
    return names


def _models(ws: Path) -> list[tuple[str, str, str]]:
    """(model type, name file, model name) for each model."""
    lines = _lines(ws / "mfsim.nam")
    begin, end = _block(lines, "MODELS")
    rows = [_tokens(ln) for ln in lines[begin + 1 : end]]
    return [(r[0].lower().removesuffix("6"), r[1], r[2]) for r in rows if r]


def _packages(ws: Path, namefile: str, mtype: str) -> list[tuple[str, str, str]]:
    """(package type, file, package name) for each package of a model. A
    package the name file gives no name is named as MODFLOW 6 names it: by
    its type, or for a type a model can have several of, `<type>-<n>`, n
    counting the packages of that type."""
    lines = _lines(ws / namefile)
    begin, end = _block(lines, "PACKAGES")
    packages, count = [], {}
    for row in (_tokens(ln) for ln in lines[begin + 1 : end]):
        if not row:
            continue
        ptype = row[0].lower().removesuffix("6")
        count[ptype] = count.get(ptype, 0) + 1
        if len(row) > 2:
            pname = row[2].lower()
        elif _multi_package(_component(mtype, ptype, _lines(ws / row[1]))):
            pname = f"{ptype}-{count[ptype]}"
        else:
            pname = ptype
        packages.append((ptype, row[1], pname))
    return packages


def _multi_package(component: str | None) -> bool:
    """Whether a model can have several packages of a DFN component."""
    if component is None:
        return False
    head = (DFN_PATH / f"{component}.dfn").read_text().splitlines()[:5]
    return any("multi-package" in line for line in head)


def _component(mtype: str, ptype: str, lines: list[str]) -> str | None:
    """The DFN component of a package file, or None if there is none."""
    words = {t.upper() for ln in lines for t in _tokens(ln)}
    variant = (
        "a" if "READASARRAYS" in words else "g" if "READARRAYGRID" in words else ""
    )
    component = f"{mtype}-{ptype}{variant}"
    return component if (DFN_PATH / f"{component}.dfn").is_file() else None


def _insert_option(lines: list[str], option: str) -> list[str]:
    span = _block(lines, "OPTIONS")
    if span is None:
        return ["BEGIN OPTIONS", f"  {option}", "END OPTIONS", ""] + lines
    return lines[: span[1]] + [f"  {option}"] + lines[span[1] :]


def _input_file(mname: str) -> str:
    return f"{mname}.export.nc"


def enable_export(ws: Path, netcdf_format: str):
    """Make each model export its package arrays to a NetCDF file in
    validate mode: the model name file names the file, and every package
    that can export its arrays gets the EXPORT_ARRAY_NETCDF option."""
    ws = Path(ws)
    for mtype, namefile, mname in _models(ws):
        packages = _packages(ws, namefile, mtype)
        if netcdf_format == "structured" and not any(p == "dis" for p, *_ in packages):
            raise ValueError(f"Structured NetCDF input requires a DIS grid: {mname}")
        nam = _lines(ws / namefile)
        if any(t.upper().startswith("NETCDF") for ln in nam for t in _tokens(ln)[:1]):
            raise ValueError(f"{namefile} already names a NetCDF file")
        option = f"{OUTPUT_OPTION[netcdf_format]} FILEOUT {_input_file(mname)}"
        (ws / namefile).write_text("\n".join(_insert_option(nam, option)) + "\n")
        for ptype, fname, _ in packages:
            lines = _lines(ws / fname)
            component = _component(mtype, ptype, lines)
            if component and "export_array_netcdf" in _dfn_names(component, "options"):
                lines = _insert_option(lines, "EXPORT_ARRAY_NETCDF")
                (ws / fname).write_text("\n".join(lines) + "\n")


def _exported(path: Path) -> set[tuple[str, str, int]]:
    """(package, tag, aux index or 0) of each array in an exported file."""
    import netCDF4

    with netCDF4.Dataset(path) as ds:
        return {
            (
                *var.modflow_input.lower().split("/")[1:3],
                int(var.modflow_iaux) if "modflow_iaux" in var.ncattrs() else 0,
            )
            for var in ds.variables.values()
            if "modflow_input" in var.ncattrs()
        }


def _aux_names(lines: list[str]) -> list[str]:
    """The package's auxiliary variable names, in order."""
    span = _block(lines, "OPTIONS")
    for line in lines[span[0] + 1 : span[1]] if span else []:
        tokens = _tokens(line)
        if tokens and tokens[0].upper() in ("AUXILIARY", "AUX"):
            return [t.lower() for t in tokens[1:]]
    return []


def _read_arrays(
    lines: list[str], block: str, names: set[str], exported: set[str]
) -> list[str]:
    """`lines` with each exported array given in a `block` block (every
    one, for PERIOD) read with NETCDF in place of its control and data
    lines. An array given another way, such as a TIMEARRAYSERIES entry,
    is kept."""
    out, inside, skip = [], False, False
    for line in lines:
        tokens = _tokens(line)
        upper = [t.upper() for t in tokens]
        if upper[:2] == ["BEGIN", block]:
            inside = True
        elif upper[:2] == ["END", block]:
            inside, skip = False, False
        elif inside and tokens and tokens[0].lower() in names:
            skip = tokens[0].lower() in exported and upper[1:] in ([], ["LAYERED"])
            if skip:
                out.append(f"  {tokens[0]} NETCDF")
                continue
        if not skip:
            out.append(line)
    return out


def use_netcdf_input(ws: Path):
    """Make each model read its exported arrays from the file
    `enable_export` had it write, instead of from its package files."""
    ws = Path(ws)
    for mtype, namefile, mname in _models(ws):
        fname = _input_file(mname)
        nam = [
            ln
            for ln in _lines(ws / namefile)
            if not (
                _tokens(ln)[:1] and _tokens(ln)[0].upper() in OUTPUT_OPTION.values()
            )
        ]
        (ws / namefile).write_text(
            "\n".join(_insert_option(nam, f"NETCDF FILEIN {fname}")) + "\n"
        )
        exported = _exported(ws / fname)
        packages = _packages(ws, namefile, mtype)
        listed = {name for ptype, _, pname in packages for name in (ptype, pname)}
        unknown = {p for p, _, _ in exported} - listed
        assert not unknown, f"{fname}: packages not in {namefile}: {sorted(unknown)}"
        n_netcdf = 0
        for ptype, pfile, pname in packages:
            original = _lines(ws / pfile)
            lines = [
                ln
                for ln in original
                if [t.upper() for t in _tokens(ln)] != ["EXPORT_ARRAY_NETCDF"]
            ]
            component = _component(mtype, ptype, lines)
            if component:
                aux = _aux_names(lines)
                tags = set()
                for p, tag, iaux in exported:
                    if p in (ptype, pname):
                        tags.add(aux[iaux - 1] if tag == "aux" and iaux else tag)
                for block in ("griddata", "period"):
                    names = _dfn_names(component, block) | set(aux)
                    lines = _read_arrays(lines, block.upper(), names, tags)
                n = sum(
                    [t.upper() for t in _tokens(ln)][1:] == ["NETCDF"] for ln in lines
                )
                assert n or not tags, f"{pfile}: no array read from {fname}"
                n_netcdf += n
            if lines != original:
                (ws / pfile).write_text("\n".join(lines) + "\n")
        assert n_netcdf, f"No package array of {namefile} read from {fname}"


# layer index arrays, which MF6 exports 1-based and flopy holds 0-based
INDEX_ARRAYS = {"irch", "ievt"}


def check_exported_arrays(ws: Path, sim) -> int:
    """
    Check the arrays each model of `sim` (a flopy simulation) exported to
    the file `enable_export` had it write hold its flopy arrays: every
    GRIDDATA array, and for array input (READASARRAYS) the first stress
    period's slice of each array given with values in that period,
    auxiliary arrays included. Arrays flopy doesn't hold as values (not
    given, or given as a time-array series) are skipped. Grid array
    (READARRAYGRID) PERIOD arrays are skipped: their slices hold values at
    boundary cells only. Returns the number of arrays checked.
    """
    import netCDF4
    import numpy as np

    n = 0
    for _, _, mname in _models(Path(ws)):
        model = sim.get_model(mname)
        with netCDF4.Dataset(Path(ws) / _input_file(mname)) as ds:
            for var in ds.variables.values():
                if "modflow_input" not in var.ncattrs():
                    continue
                _, pname, tag = var.modflow_input.lower().split("/")
                iaux = int(var.modflow_iaux) if "modflow_iaux" in var.ncattrs() else 0
                package = model.get_package(pname)
                assert package is not None, f"{mname}: no package {pname} in flopy"
                data = getattr(package, tag, None)
                if data is None:
                    continue
                values = var[:]
                if "time" in var.dimensions:
                    readasarrays = getattr(package, "readasarrays", None)
                    if readasarrays is None or not readasarrays.get_data():
                        continue
                    try:
                        expected = data.get_data(0)
                    except Exception:
                        continue
                    values = values[0]
                    if iaux and expected is not None:
                        expected = np.asarray(expected)[iaux - 1]
                    if tag in INDEX_ARRAYS and expected is not None:
                        expected = np.asarray(expected) + 1
                else:
                    expected = data.array
                if expected is None:
                    continue
                expected = np.asarray(expected)
                if "layer" in var.ncattrs():
                    expected = expected[int(var.layer) - 1]
                assert np.allclose(np.ravel(values), np.ravel(expected)), (
                    f"{mname}: exported {var.name} differs from {pname}.{tag}"
                )
                n += 1
    return n


# model type: (dependent variable, its output control record word)
DEPENDENT = {
    "gwf": ("head", "HEAD"),
    "gwt": ("concentration", "CONCENTRATION"),
    "gwe": ("temperature", "TEMPERATURE"),
}


def _output_file(mname: str) -> str:
    return f"{mname}.nc"


def enable_output(ws: Path, netcdf_format: str):
    """Make each model also write its dependent variable to a NetCDF output
    file, in the format its input was exported in."""
    ws = Path(ws)
    for _, namefile, mname in _models(ws):
        option = f"{OUTPUT_OPTION[netcdf_format]} FILEOUT {_output_file(mname)}"
        nam = _insert_option(_lines(ws / namefile), option)
        (ws / namefile).write_text("\n".join(nam) + "\n")


def check_outputs(ws: Path) -> int:
    """Check each model's NetCDF output file `enable_output` had it write
    against the binary output file its output control saves, if any.
    Returns the number of models checked."""
    from netcdf_output import check_netcdf_output

    ws = Path(ws)
    n = 0
    for mtype, namefile, mname in _models(ws):
        if mtype not in DEPENDENT:
            continue
        text, word = DEPENDENT[mtype]
        oc = [f for p, f, _ in _packages(ws, namefile, mtype) if p == "oc"]
        binary = None
        for line in _lines(ws / oc[0]) if oc else []:
            tokens = _tokens(line)
            if [t.upper() for t in tokens[:2]] == [word, "FILEOUT"]:
                binary = tokens[2].strip("'\"")
        if binary is None:
            continue
        check_netcdf_output(ws / _output_file(mname), ws / binary, text)
        n += 1
    return n

import argparse
import textwrap
from collections import Counter
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from pprint import pprint
from typing import Optional

from filters import Filters
from jinja2 import Environment, FileSystemLoader
from modflow_devtools.dfn import Dfn

MF6_LENVARNAME = 16
F90_LINELEN = 82
PROJ_ROOT_PATH = Path(__file__).parents[3]
TEMPLATES_PATH = Path(__file__).parents[1] / "templates"
DEFAULT_DFNS_PATH = Path(__file__).parents[1] / "dfns.txt"
DFN_PATH = PROJ_ROOT_PATH / "doc" / "mf6io" / "mf6ivar" / "dfn"
SRC_PATH = PROJ_ROOT_PATH / "src"
IDM_PATH = SRC_PATH / "Idm"
# overhead chars for the two write-statement forms used by spec_line
_WRITE_STD = "    write(output_unit, '(a)') "
_WRITE_NO_ADV = "    write(output_unit, '(a)', advance='no') "
_MAX_STD = 132 - len(_WRITE_STD) - 2
_MAX_NO_ADV = 132 - len(_WRITE_NO_ADV) - 2


def _spec_line(s: str) -> str:
    """Jinja filter: emit Fortran write statement(s) for one DFN text line.

    Short lines produce a single write; long lines are split into
    advance='no' chunks so generated source lines stay within 132 chars.
    Single-quotes in content are doubled to form valid Fortran literals.
    Chunks never split a '' pair: if a slice ends with an odd number of
    trailing quotes, it is shortened by one character so the split always
    falls between pairs.
    """
    escaped = s.replace("'", "''")
    if len(escaped) <= _MAX_STD:
        return f"{_WRITE_STD}'{escaped}'"
    chunks = []
    while escaped:
        chunk = escaped[:_MAX_NO_ADV]
        n_trailing = len(chunk) - len(chunk.rstrip("'"))
        if n_trailing % 2 == 1:
            chunk = escaped[: _MAX_NO_ADV - 1]
        chunks.append(chunk)
        escaped = escaped[len(chunk) :]
    lines = [f"{_WRITE_NO_ADV}'{c}'" for c in chunks]
    lines.append("    write(output_unit, '(a)') ''")
    return "\n".join(lines)


_BASE_TYPE_MAP = {
    "double precision": "DOUBLE",
    "integer": "INTEGER",
    "keyword": "KEYWORD",
    "string": "STRING",
}


def _normalize_type(t_raw: str, shape_str: str, ndim: int, aggregate: bool) -> str:
    """Map a raw DFN type string to its IDM Fortran representation."""
    if aggregate:
        return t_raw.upper()
    t_lower = t_raw.lower()
    if t_lower in _BASE_TYPE_MAP:
        t = _BASE_TYPE_MAP[t_lower]
        if shape_str and t in ("DOUBLE", "INTEGER"):
            t = f"{t}{ndim}D"
        return t
    return t_raw.upper()


def _wrap_f90_content(raw: str, max_width: int = 60) -> str:
    """Wrap a string for use in a Fortran string literal.

    Returns the raw text split into newline-separated segments, each at
    most ``max_width`` characters wide.  The ``| value`` Jinja filter
    detects the embedded newlines and emits ``// &`` string concatenation
    so that the generated Fortran source lines stay within the project
    line-length limit.  Single-segment strings are returned as-is.
    """
    llist = textwrap.wrap(raw, max_width)
    if not llist:
        return ""
    if len(llist) == 1:
        return llist[0]
    # Preserve the inter-word space at each line break: the space between
    # the last word of one segment and the first word of the next is dropped
    # by textwrap.wrap, so add it back as a trailing space on each non-last
    # segment.  The | value filter will join quoted segments with //.
    parts = [ln + " " for ln in llist[:-1]] + [llist[-1]]
    return "\n".join(parts)


@dataclass
class Param:
    """Represents a single input parameter definition from a DFN file."""

    component: str
    subcomponent: str
    block: str
    tag: str
    fortran_var: str
    type: str
    shape: str
    longname: str
    required: bool
    developmode: bool
    in_record: bool
    preserve_case: bool
    layered: bool
    timeseries: bool
    aggregate: bool = False
    block_variable: bool = False
    block_qualified: bool = False

    @property
    def varname(self) -> str:
        """Full Fortran variable name for this parameter definition."""
        if self.block_qualified:
            return (
                f"{self.component.lower()}"
                f"{self.subcomponent.lower()}"
                f"_{self.block.lower()}"
                f"_{self.fortran_var.lower()}"
            )
        return (
            f"{self.component.lower()}"
            f"{self.subcomponent.lower()}"
            f"_{self.fortran_var.lower()}"
        )

    @property
    def found_name(self) -> str:
        """Field name used in the ParamFoundType derived type."""
        if self.block_qualified:
            return f"{self.block.lower()}_{self.fortran_var.lower()}"
        return self.fortran_var.lower()


@dataclass
class Block:
    """Represents an input block from a DFN file."""

    name: str
    required: bool
    aggregate: bool
    block_var: bool


@dataclass
class DfnFile:
    """Parsed representation of a DFN file."""

    component: str
    subcomponent: str
    multi_package: bool
    subpackages: list
    params: list  # all params (aggregate + non-aggregate), excluding block_variable
    blocks: list

    @property
    def param_definitions(self) -> list:
        """Non-aggregate params for param_definitions array."""
        return [p for p in self.params if not p.aggregate]

    @property
    def aggregate_definitions(self) -> list:
        """Aggregate params for aggregate_definitions array."""
        return [p for p in self.params if p.aggregate]


def parse_dfn(dfnfspec: Path, common: Optional[dict] = None) -> DfnFile:
    """Parse a DFN file into a DfnFile object."""
    component, subcomponent = dfnfspec.stem.upper().split("-")

    # Pre-scan for multi_package and mf6 subpackages.
    # _load_v1_flat only captures "# flopy ..." comments, not "# mf6 subpackage".
    multi_package = False
    subpackages = []
    for line in dfnfspec.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped.startswith("#"):
            continue
        if "flopy multi-package" in stripped:
            multi_package = True
        elif "mf6 subpackage" in stripped:
            sp = stripped.replace("# mf6 subpackage", "").strip().upper()
            subpackages.append(sp.ljust(16))

    if not subpackages:
        subpackages = [" " * 16]

    # Parse variable entries using modflow_devtools
    with dfnfspec.open(encoding="utf-8") as f:
        flat, _ = Dfn._load_v1_flat(f, common=common)

    # Track blocks in DFN order
    block_names_ordered = []
    block_data = {}  # blockname -> tracking dict

    params = []

    for vd in flat.values(multi=True):
        blockname = vd.get("block", "")
        if not blockname:
            continue

        blockname_upper = blockname.upper()

        if blockname_upper not in block_data:
            block_names_ordered.append(blockname_upper)
            block_data[blockname_upper] = {
                "required_l": [],
                "has_block_var": False,
                "is_aggregate": False,
                "aggregate_required": False,
            }

        is_block_variable = vd.get("block_variable", "").lower() == "true"
        if is_block_variable:
            block_data[blockname_upper]["has_block_var"] = True
            continue

        vn = vd["name"].upper()
        mf6vn = vd["mf6internal"].upper() if "mf6internal" in vd else vn
        # sanitize: hyphens are not valid in Fortran identifiers
        mf6vn = mf6vn.replace("-", "_")

        t_raw = vd.get("type", "")
        aggregate_t = t_raw.lower().startswith("recarray") or t_raw.lower().startswith(
            "keystring"
        )

        # Shape processing
        shape = vd.get("shape", "")
        if component.upper() == "EXG" and vn in ("CELLIDM1", "CELLIDM2"):
            shape = "(ncelldim)"
        shape = shape.replace("(", "").replace(")", "").replace(",", "").upper()
        if mf6vn == "AUXVAR":
            if shape == "NCOL*NROW; NCPL":
                shape = "NAUX NCPL"
            elif shape == "NODES":
                shape = "NAUX NODES"
        elif shape == "NCOL*NROW; NCPL":
            shape = "NCPL"
        shapelist = shape.strip().split() if shape.strip() else []
        ndim = len(shapelist)
        shape_str = " ".join(shapelist)

        t = _normalize_type(t_raw, shape_str, ndim, aggregate_t)
        if len(t) > 60:
            t = _wrap_f90_content(t, max_width=60)

        # Longname wrapping
        longname = ""
        if vd.get("longname"):
            raw = vd["longname"].replace("'", "")
            longname = _wrap_f90_content(raw, max_width=60)

        required = vd.get("optional", "").lower() != "true"
        developmode = vd.get("developmode", "").lower() == "true"
        in_record = vd.get("in_record", "").lower() == "true"
        preserve_case = vd.get("preserve_case", "").lower() == "true"
        layered = vd.get("layered", "").lower() == "true"
        timeseries = vd.get("time_series", "").lower() == "true"

        param = Param(
            component=component,
            subcomponent=subcomponent,
            block=blockname_upper,
            tag=vn,
            fortran_var=mf6vn,
            type=t,
            shape=shape_str,
            longname=longname,
            required=required,
            developmode=developmode,
            in_record=in_record,
            preserve_case=preserve_case,
            layered=layered,
            timeseries=timeseries,
            aggregate=aggregate_t,
        )
        params.append(param)

        # Update block tracking
        if aggregate_t:
            block_data[blockname_upper]["is_aggregate"] = True
            block_data[blockname_upper]["aggregate_required"] = required
        elif not in_record:
            block_data[blockname_upper]["required_l"].append(required)

    # Prepend OPTIONS block if not already present
    if block_names_ordered and "OPTIONS" not in block_names_ordered:
        block_names_ordered.insert(0, "OPTIONS")
        block_data["OPTIONS"] = {
            "required_l": [],
            "has_block_var": False,
            "is_aggregate": False,
            "aggregate_required": False,
        }

    # --- Keystring member expansion ---
    # For blocks containing a keystring param, only the keystring placeholder
    # itself (e.g. MAWSETTING, TVKSETTING) is excluded from param_definitions.
    # All keystring member params are retained so the Fortran IDM framework
    # can dispatch them.
    keystring_blocks: set = set()  # blocks that contain a keystring param
    for p in params:
        if p.type.upper().startswith("KEYSTRING"):
            keystring_blocks.add(p.block)

    if keystring_blocks:
        # Collect the keystring placeholder tag(s) to exclude per block
        excl: dict = {}
        for p in params:
            if p.block not in keystring_blocks:
                continue
            if p.type.upper().startswith("KEYSTRING"):
                excl.setdefault(p.block, set()).add(p.tag)

        params = [
            p for p in params if p.aggregate or p.tag not in excl.get(p.block, set())
        ]

    # --- Detect fortran_var collisions and mark block-qualified names ---
    # When the same fortran_var name appears in more than one block, prefix
    # the generated Fortran constant name and ParamFoundType field with the
    # block name to avoid duplicate Fortran identifiers.
    var_counts = Counter(p.fortran_var.upper() for p in params)
    for p in params:
        if var_counts[p.fortran_var.upper()] > 1:
            p.block_qualified = True

    # Build Block objects in DFN order
    blocks = []
    for blockname_upper in block_names_ordered:
        bdata = block_data[blockname_upper]
        if bdata["is_aggregate"]:
            block_required = bdata["aggregate_required"]
        else:
            block_required = any(bdata["required_l"])

        blocks.append(
            Block(
                name=blockname_upper,
                required=block_required,
                aggregate=bdata["is_aggregate"],
                block_var=bdata["has_block_var"],
            )
        )

    return DfnFile(
        component=component,
        subcomponent=subcomponent,
        multi_package=multi_package,
        subpackages=subpackages,
        params=params,
        blocks=blocks,
    )


def _get_template_env() -> Environment:
    template_loader = FileSystemLoader(str(TEMPLATES_PATH))
    template_env = Environment(
        loader=template_loader,
        trim_blocks=True,
        lstrip_blocks=True,
        keep_trailing_newline=True,
    )
    template_env.filters["value"] = Filters.value
    template_env.filters["spec_line"] = _spec_line
    return template_env


def make_spec(outpath: PathLike | None = None, verbose: bool = False) -> None:
    """Generate DfnSpec.f90, writing to outpath or stdout if outpath is None."""
    import sys
    import warnings

    template_env = _get_template_env()
    template = template_env.get_template("DfnSpec.f90.jinja")

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*modflow_devtools.dfns.*experimental.*"
        )
        from modflow_devtools.dfns import DfnSpec

        spec = DfnSpec.load(DFN_PATH)

    lines = spec.dumps().splitlines()
    content = template.render(lines=lines)

    if outpath is None:
        sys.stdout.write(content)
    else:
        if verbose:
            print(f"  writing {outpath}")
        with open(outpath, "w", newline="\n") as f:
            f.write(content)


def make_targets(dfn: DfnFile, outdir: PathLike, verbose: bool = False):
    """Generate the Fortran IDM file for a single parsed DfnFile."""
    outdir = Path(outdir)
    template_env = _get_template_env()
    component_idm_template = template_env.get_template("Componentidm.f90.jinja")
    ofname = f"{dfn.component.lower()}-{dfn.subcomponent.lower()}idm.f90"
    ofspec = outdir / ofname
    if verbose:
        print(f"  writing {ofspec}")
    with open(ofspec, "w", newline="\n") as f:
        f.write(component_idm_template.render(dfn=dfn))


def make_all(
    dfn_paths: list,
    outdir: PathLike,
    verbose: bool = False,
):
    """Generate all Fortran IDM and selector source files from a list of DFN paths."""
    outdir = Path(outdir)
    selector_dir = outdir / "selector"
    selector_dir.mkdir(parents=True, exist_ok=True)

    template_env = _get_template_env()
    component_idm_template = template_env.get_template("Componentidm.f90.jinja")
    component_selector_template = template_env.get_template(
        "IdmComponentDfnSelector.f90.jinja"
    )
    selector_template = template_env.get_template("IdmDfnSelector.f90.jinja")

    # Load common variables for description substitution
    common_path = DFN_PATH / "common.dfn"
    common = None
    if common_path.is_file():
        if verbose:
            print(f"  loading {common_path}")
        with common_path.open(encoding="utf-8") as f:
            common, _ = Dfn._load_v1_flat(f)

    # Parse all DFN files
    dfn_files = []
    for path in dfn_paths:
        if verbose:
            print(f"  parsing {path}")
        dfn_files.append(parse_dfn(path, common=common))

    # Write component IDM files
    for dfn in dfn_files:
        ofname = f"{dfn.component.lower()}-{dfn.subcomponent.lower()}idm.f90"
        ofspec = outdir / ofname
        if verbose:
            print(f"  writing {ofspec}")
        with open(ofspec, "w", newline="\n") as f:
            f.write(component_idm_template.render(dfn=dfn))

    # Group DFN files by component, preserving within-component order
    components = {}
    for dfn in dfn_files:
        if dfn.component not in components:
            components[dfn.component] = []
        components[dfn.component].append(dfn)

    # Write component selector files
    for component, packages in components.items():
        ofname = f"Idm{component.title()}DfnSelector.f90"
        ofspec = selector_dir / ofname
        if verbose:
            print(f"  writing {ofspec}")
        with open(ofspec, "w", newline="\n") as f:
            f.write(
                component_selector_template.render(
                    component=component,
                    packages=packages,
                )
            )

    # Write master selector file
    ofspec = selector_dir / "IdmDfnSelector.f90"
    if verbose:
        print(f"  writing {ofspec}")
    with open(ofspec, "w", newline="\n") as f:
        f.write(selector_template.render(components=list(components.keys())))


def _expand_dfns(dfns_arg) -> list:
    """
    Expand DFN file paths, a dfns.txt listing,
    or a directory to a list of DFN Paths.
    """
    if isinstance(dfns_arg, list):
        # Strip whitespace (including \r from Windows line endings) from each path
        paths = [Path(str(p).strip()) for p in dfns_arg]
    elif isinstance(dfns_arg, (str, Path)):
        paths = [Path(str(dfns_arg).strip())]
    else:
        raise TypeError(f"Unexpected type: {type(dfns_arg)}")

    result = []
    for path in paths:
        if path.is_dir():
            result.extend(sorted(path.glob("*.dfn")))
        elif path.is_file() and path.suffix == ".txt":
            # Read list of filenames from text file
            for fname in path.read_text(encoding="utf-8").splitlines():
                fname = fname.strip()
                if fname and not fname.startswith("#"):
                    p = DFN_PATH / fname
                    if p.is_file():
                        result.append(p)
        else:
            # Check if it's a simple filename that needs DFN_PATH prepended
            if len(path.parts) == 1:
                path = DFN_PATH / path
            if path.exists() and path.is_file():
                result.append(path)

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Convert DFN files to Fortran source files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Generate Fortran source code from DFN files. This script
            converts definition (DFN) files to Fortran source files,
            each representing a parameter set for a particular input
            definition. Fortran files generated by this tool provide
            support for simulations, models or packages described by
            the given DFN files. Each DFN file is transformed into a
            corresponding Fortran file with "idm" and the same stem:
            e.g. gwf-ic.dfn becomes gwf-icidm.f90.
            """
        ),
    )
    parser.add_argument(
        "dfn",
        nargs="*",
        default=DEFAULT_DFNS_PATH,
        help="Path(s) to DFN files, a dfns.txt listing, or a directory of DFN files",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=False,
        default=IDM_PATH,
        help="The directory to write Fortran source files",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        required=False,
        default=False,
        help="Whether to show verbose output",
    )
    parser.add_argument(
        "--spec",
        metavar="OUTPUT",
        nargs="?",
        const="-",
        default=None,
        help="Generate DfnSpec.f90 only, writing to OUTPUT ('-' or omitted: stdout)",
    )
    args = parser.parse_args()
    verbose = args.verbose

    if args.spec is not None:
        outpath = None if args.spec == "-" else Path(args.spec)
        make_spec(outpath=outpath, verbose=verbose)
        raise SystemExit(0)

    dfn_arg = args.dfn if args.dfn else DEFAULT_DFNS_PATH
    dfn_paths = _expand_dfns(dfn_arg)
    outdir = Path(args.outdir) if args.outdir else Path.cwd()

    if verbose:
        print("Generating Fortran source files from DFNs:")
        pprint(dfn_paths)

    make_all(dfn_paths, outdir, verbose=verbose)

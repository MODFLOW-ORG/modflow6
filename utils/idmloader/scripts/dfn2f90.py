import argparse
import textwrap
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from pprint import pprint

from filters import Filters
from jinja2 import Environment, FileSystemLoader

MF6_LENVARNAME = 16
F90_LINELEN = 82
PROJ_ROOT_PATH = Path(__file__).parents[3]
TEMPLATES_PATH = Path(__file__).parents[1] / "templates"
DEFAULT_DFNS_PATH = Path(__file__).parents[1] / "dfns.txt"
DFN_PATH = PROJ_ROOT_PATH / "doc" / "mf6io" / "mf6ivar" / "dfn"
SRC_PATH = PROJ_ROOT_PATH / "src"
IDM_PATH = SRC_PATH / "Idm"


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

    @property
    def varname(self) -> str:
        """Full Fortran variable name for this parameter definition."""
        return (
            f"{self.component.lower()}"
            f"{self.subcomponent.lower()}"
            f"_{self.fortran_var.lower()}"
        )


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


def parse_dfn(dfnfspec: Path) -> DfnFile:
    """Parse a DFN file into a DfnFile object."""
    component, subcomponent = dfnfspec.stem.upper().split("-")
    multi_package = False
    subpackages = []

    lines = dfnfspec.read_text(encoding="utf-8").splitlines()

    # Parse raw variable entries in DFN file order
    var_entries = []  # list of ((name, block), dict) to preserve order
    vd = {}

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if vd:
                name = vd.get("name", "")
                block = vd.get("block", "")
                var_entries.append(((name, block), vd))
            vd = {}
            continue

        if stripped.startswith("#"):
            if "flopy multi-package" in stripped:
                multi_package = True
            elif "mf6 subpackage" in stripped:
                sp = stripped.replace("# mf6 subpackage", "").strip().upper()
                subpackages.append(sp.ljust(16))
            continue

        parts = stripped.split(None, 1)
        if parts:
            k = parts[0]
            v = parts[1] if len(parts) > 1 else ""
            vd[k] = v

    if vd:
        name = vd.get("name", "")
        block = vd.get("block", "")
        var_entries.append(((name, block), vd))

    if not subpackages:
        subpackages = [" " * 16]

    # Track blocks in DFN order
    block_names_ordered = []
    block_data = {}  # blockname -> tracking dict

    params = []

    for (name, blockname), vd in var_entries:
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

        is_block_variable = vd.get("block_variable", "").upper() == "TRUE"
        if is_block_variable:
            block_data[blockname_upper]["has_block_var"] = True
            continue

        vn = name.upper()
        mf6vn = vd["mf6internal"].upper() if "mf6internal" in vd else vn

        t = vd.get("type", "").upper()
        aggregate_t = t.startswith("RECARRAY")

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

        # Type normalization (non-aggregate only)
        if not aggregate_t:
            if t == "DOUBLE PRECISION":
                t = "DOUBLE"
            if shape_str and t in ("DOUBLE", "INTEGER"):
                t = f"{t}{ndim}D"

        # Longname wrapping
        longname = ""
        if vd.get("longname"):
            raw = vd["longname"].replace("'", "")
            llist = textwrap.wrap(raw, 70)
            if len(llist) == 1:
                longname = llist[0]
            elif len(llist) > 1:
                longname = f"{llist[0]}&\n"
                for ln in llist[1:-1]:
                    longname += f"     & {ln}&\n"
                longname += f"     & {llist[-1]}"

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
        line_statement_prefix="_",
        keep_trailing_newline=True,
    )
    template_env.filters["value"] = Filters.value
    return template_env


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

    # Parse all DFN files
    dfn_files = []
    for path in dfn_paths:
        if verbose:
            print(f"  parsing {path}")
        dfn_files.append(parse_dfn(path))

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
        paths = [Path(p) for p in dfns_arg]
    elif isinstance(dfns_arg, (str, Path)):
        paths = [Path(dfns_arg)]
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
        elif path.is_file():
            if len(path.parts) == 1:
                path = DFN_PATH / path
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
    args = parser.parse_args()
    dfn_arg = args.dfn if args.dfn else DEFAULT_DFNS_PATH
    dfn_paths = _expand_dfns(dfn_arg)
    outdir = Path(args.outdir) if args.outdir else Path.cwd()
    verbose = args.verbose

    if verbose:
        print("Generating Fortran source files from DFNs:")
        pprint(dfn_paths)

    make_all(dfn_paths, outdir, verbose=verbose)

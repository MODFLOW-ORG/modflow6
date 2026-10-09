import os
import shutil
from collections.abc import Callable, Iterable
from itertools import repeat
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from traceback import format_exc
from warnings import warn

import flopy
from compare import (
    Comparison,
    adjust_htol,
    detect_comparison,
    get_comparison_files,
    get_mf6_files,
    get_namefiles,
    get_rclose,
)
from flopy.mbase import BaseModel
from flopy.mf6 import MFSimulation
from flopy.utils.compare import compare_cell_budget, compare_heads
from modflow_devtools.misc import get_ostag, is_in_ci

DNODATA = 3.0e30
EXTTEXT = {
    "hds": "head",
    "hed": "head",
    "bhd": "head",
    "ucn": "concentration",
    "cbc": "cell-by-cell",
}
HDS_EXT = (
    "hds",
    "hed",
    "bhd",
    "ahd",
    "bin",
)
CBC_EXT = (
    "cbc",
    "bud",
)
# model output, not copied to a re-written simulation's workspace
OUTPUT_EXT = {
    ".lst",
    ".hds",
    ".hed",
    ".bhd",
    ".ahd",
    ".ucn",
    ".cbc",
    ".bud",
    ".grb",
    ".csv",
    ".out",
}
FP4_COMPARISONS = (Comparison.FP4_STRUCTURED, Comparison.FP4_LAYERED)


def api_return(success, model_ws) -> tuple[bool, list[str]]:
    """
    parse libmf6 stdout shared object file
    """
    fpth = os.path.join(model_ws, "mfsim.stdout")
    return success, open(fpth).readlines()


def get_mfsim_lst_tail(path: os.PathLike, lines=100) -> str:
    """Get the tail of the mfsim.lst listing file"""
    msg = ""
    _lines = open(path).read().splitlines()
    msg = "\n" + 79 * "-" + "\n"
    i0 = -lines if len(_lines) > lines else 0
    for line in _lines[i0:]:
        if len(line) > 0:
            msg += f"{line}\n"
    msg += 79 * "-" + "\n\n"
    return msg


def get_workspace(sim_or_model) -> Path:
    if isinstance(sim_or_model, MFSimulation):
        return sim_or_model.sim_path
    elif isinstance(sim_or_model, BaseModel):
        return Path(sim_or_model.model_ws)
    else:
        raise ValueError(f"Unsupported model type: {type(sim_or_model)}")


def run_parallel(workspace, target, ncpus) -> tuple[bool, list[str]]:
    if not is_in_ci() and get_ostag() in ["mac"]:
        oversubscribed = ["--hostfile", "localhost"]
        with open(f"{workspace}/localhost", "w") as f:
            f.write(f"localhost slots={ncpus}\n")
    else:
        oversubscribed = ["--oversubscribe"]

    normal_msg = "normal termination"
    success = False
    nr_success = 0
    buff = []

    # parallel commands
    if get_ostag() in ["win64"]:
        mpiexec_cmd = ["mpiexec", "-np", str(ncpus), target, "-p"]
    else:
        mpiexec_cmd = ["mpiexec"] + oversubscribed + ["-np", str(ncpus), target, "-p"]

    proc = Popen(mpiexec_cmd, stdout=PIPE, stderr=STDOUT, cwd=workspace)

    while True:
        line = proc.stdout.readline().decode("utf-8")
        if line == "" and proc.poll() is not None:
            break
        if line:
            # success is when the success message appears
            # in every process of the parallel simulation
            if normal_msg in line.lower():
                nr_success += 1
                if nr_success == ncpus:
                    success = True
            line = line.rstrip("\r\n")
            print(line)
            buff.append(line)
        else:
            break

    return success, buff


def write_input(*sims, overwrite: bool = True, verbose: bool = True):
    """
    Write input files for `flopy.mf6.MFSimulation` or `flopy.mbase.BaseModel`.

    Parameters
    ----------

    sims : arbitrary list
        Simulations or models
    verbose : bool, optional
        whether to show verbose output
    """

    if sims is None:
        warn("No simulations or models!")
        return

    # write input files for each model or simulation
    for sim in sims:
        if sim is None:
            continue

        if isinstance(sim, flopy.mf6.MFSimulation):
            workspace = Path(sim.sim_path)
            if any(workspace.glob("*")) and not overwrite:
                warn("Workspace is not empty, not writing input files")
                return
            if verbose:
                print(f"Writing mf6 simulation '{sim.name}' to: {sim.sim_path}")
            sim.write_simulation()
        elif isinstance(sim, flopy.mbase.BaseModel):
            workspace = Path(sim.model_ws)
            if any(workspace.glob("*")) and not overwrite:
                warn("Workspace is not empty, not writing input files")
                return
            if verbose:
                print(f"Writing {type(sim)} model '{sim.name}' to: {sim.model_ws}")
            sim.write_input()
        else:
            raise ValueError(f"Unsupported simulation/model type: {type(sim)}")


def dependent_variable_text(fpth, default="concentration") -> str:
    """Text of the dependent variable records in a binary output file
    (e.g. GWT concentration or GWE temperature), for flopy's readers."""
    for text in (default, "concentration", "temperature", "head"):
        try:
            flopy.utils.HeadFile(fpth, precision="double", text=text)
            return text
        except Exception:
            continue
    return default


def require_flopy4():
    """Skip a flopy4 comparison without flopy4, except in CI, where it fails."""
    import pytest

    try:
        import flopy4  # noqa: F401
    except ImportError:
        if is_in_ci():
            raise
        pytest.skip("flopy4 not installed")


def write_fp4_netcdf(
    src: os.PathLike,
    dst: os.PathLike,
    netcdf_format: str,
    sim_dirs: Iterable[os.PathLike] = (".",),
):
    """
    Re-write the simulations in `src` into `dst` with flopy4, with each
    model's package array input in a NetCDF file.

    Parameters
    ----------
    src : path-like
        Workspace of the test to re-write.
    dst : path-like
        Workspace to write to (replaced if it exists).
    netcdf_format : str
        "structured" (DIS grids only) or "layered" (UGRID layered mesh).
    sim_dirs : iterable of path-like
        Simulation workspaces to re-write, relative to `src`.
    """
    from flopy4.mf6 import NetCDFFormat, Simulation
    from flopy4.mf6.netcdf import NetCDFModel
    from flopy4.mf6.write_context import WriteContext

    src, dst = Path(src), Path(dst)
    sim_dirs = [Path(d) for d in sim_dirs]
    # copy all but the simulations' outputs and other comparisons; this keeps
    # files flopy4 doesn't load (e.g. OBS, TS)
    skip_dirs = {c.value for c in Comparison} - {
        d.parts[0] for d in sim_dirs if d.parts
    }

    # output files each simulation names (e.g. a .bin head file)
    named_outputs = {
        d: {Path(f).name for f in get_mf6_files(src / d / "mfsim.nam")[1]}
        for d in sim_dirs
    }

    def ignore(dirpath, names):
        rel = Path(dirpath).relative_to(src)
        ignored = {n for n in names if rel == Path(".") and n in skip_dirs}
        if rel in sim_dirs:
            ignored |= {
                n
                for n in names
                if (Path(n).suffix.lower() in OUTPUT_EXT or n in named_outputs[rel])
                and (Path(dirpath) / n).is_file()
            }
        return ignored

    if dst.is_dir():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=ignore)

    fmt = NetCDFFormat(netcdf_format)
    for sim_dir in sim_dirs:
        sim = Simulation.load(dst / sim_dir / "mfsim.nam")
        for name, model in sim.models.items():
            if fmt == NetCDFFormat.STRUCTURED and type(model.dis).__name__ != "Dis":
                raise ValueError(f"Structured NetCDF input requires a DIS grid: {name}")
            # the name file refers to it relative to the workspace
            nc_name = f"{name}.input.nc"
            model.netcdf_input_file = Path(nc_name)
            NetCDFModel.from_model(model, netcdf_format=fmt).to_netcdf(
                dst / sim_dir / nc_name
            )
        with WriteContext(use_netcdf=True):
            sim.write()


DFN_PATH = Path(__file__).parents[1] / "doc" / "mf6io" / "mf6ivar" / "dfn"


def _rows(path: Path) -> list[list[str]]:
    """Non-comment rows of an input file, split into tokens."""
    rows = [ln.split() for ln in path.read_text(errors="ignore").splitlines()]
    return [r for r in rows if r and not r[0].startswith(("#", "!"))]


def _block(rows: list[list[str]], name: str) -> list[list[str]]:
    """The rows of the first block named `name`."""
    upper = [[t.upper() for t in r] for r in rows]
    begin = upper.index(["BEGIN", name])
    end = upper.index(["END", name], begin)
    return rows[begin + 1 : end]


def netcdf_arrays(component: str) -> set[str]:
    """Names of the arrays of a DFN component that MF6 can read from NetCDF."""
    names, name = set(), None
    for row in (DFN_PATH / f"{component}.dfn").read_text().splitlines():
        tokens = row.split()
        if tokens[:1] == ["name"] and len(tokens) > 1:
            name = tokens[1]
        elif tokens == ["netcdf", "true"] and name:
            names.add(name)
    return names


def package_arrays(path: Path) -> list[tuple[str, bool]]:
    """(name, read from NetCDF) for each array in a package file."""
    rows = _rows(path)
    arrays = []
    for i, row in enumerate(rows):
        upper = [t.upper() for t in row]
        if (
            len(upper) > 1
            and upper[-1] == "NETCDF"
            and upper[0] not in ("BEGIN", "END")
        ):
            arrays.append((row[0].lower(), True))
        elif (
            upper[1:] in ([], ["LAYERED"])
            and i + 1 < len(rows)
            and rows[i + 1][0].upper() in ("CONSTANT", "INTERNAL", "OPEN/CLOSE")
        ):
            arrays.append((row[0].lower(), False))
    return arrays


def check_fp4_netcdf_input(workspace: os.PathLike):
    """
    Check a simulation written by `write_fp4_netcdf` reads every array MF6
    can read from NetCDF (per the DFNs) with NETCDF, from a variable in its
    model's NetCDF4 input file.
    """
    import netCDF4

    workspace = Path(workspace)
    models = _block(_rows(workspace / "mfsim.nam"), "MODELS")
    assert models, f"No models in {workspace / 'mfsim.nam'}"
    for mtype, namefile, mname in (r[:3] for r in models):
        mtype = mtype.lower().removesuffix("6")
        rows = _rows(workspace / namefile)
        nc = [r[2] for r in rows if [t.upper() for t in r[:2]] == ["NETCDF", "FILEIN"]]
        assert nc, f"No NETCDF FILEIN in {namefile}"
        with netCDF4.Dataset(workspace / nc[0]) as ds:
            assert ds.data_model == "NETCDF4", f"{nc[0]}: {ds.data_model}"
            assert "modflow_model" in ds.ncattrs(), f"{nc[0]}: no modflow_model"
            inputs = {
                ds[v].modflow_input.lower()
                for v in ds.variables
                if "modflow_input" in ds[v].ncattrs()
            }
        n_netcdf = 0
        for ftype, fname, *pname in _block(rows, "PACKAGES"):
            ptype = ftype.lower().removesuffix("6")
            # variables are named by package type, or by package name for
            # types a model can have several of
            pnames = {ptype, pname[0].lower() if pname else ptype}
            pkg_rows = _rows(workspace / fname)
            options = [t.upper() for r in pkg_rows for t in r]
            variant = (
                "a"
                if "READASARRAYS" in options
                else "g"
                if "READARRAYGRID" in options
                else ""
            )
            component = f"{mtype}-{ptype}{variant}"
            if not (DFN_PATH / f"{component}.dfn").is_file():
                continue
            capable = netcdf_arrays(component)
            aux = {
                t.lower()
                for r in pkg_rows
                if r[0].upper() in ("AUXILIARY", "AUX")
                for t in r[1:]
            }
            for name, from_netcdf in package_arrays(workspace / fname):
                tag = "aux" if name in aux else name
                if from_netcdf:
                    n_netcdf += 1
                    keys = {f"{mname}/{p}/{tag}".lower() for p in pnames}
                    assert keys & inputs, f"{fname}: {name} has no variable in {nc[0]}"
                else:
                    assert tag not in capable, f"{fname}: {name} not read from NetCDF"
        assert n_netcdf, f"No package array of {namefile} read from NetCDF"


class TestFramework:
    """
    Defines a MODFLOW 6 end-to-end integration test. Configurable
    hooks are available to build models and evaluate/plot results.
    Regression testing is also supported, allowing for comparisons
    with previous versions of MODFLOW 6 or other MODFLOW programs.

    There are three main hooks:

        - `build`: function to build the simulation(s) or model(s)
        - `check`: function to evaluate the results of the model(s)
        - `plot`: function to create plots of the results

    There are several supported comparison scenarios: typically,
    the test framework compares the output of the MF6 under test
    to the latest MF6 release. The "original" regression testing
    approach is also still supported, in which MF6 is compared to
    another program, e.g. MODFLOW-2005, MODFLOW-NWT, MODFLOW-USG.

    Parameters
    ----------
    name : str
        The test name
    workspace : pathlike
        The test workspace
    targets : dict
        Binary targets to test against. Development binaries are
        required, downloads/rebuilt binaries are optional (if not
        found, comparisons and regression tests will be skipped).
        Dictionary maps target names to paths. The test framework
        will refuse to run a program if it is not a known target.
    api_func: callable, optional
        User defined function to invoke the MODFLOW API, accepting
        the MF6 library path and the test workspace as parameters.
    build : callable, optional
        User defined function returning one or more simulations/models.
        Takes `self` as input. This is the place to build simulations.
        If no build function is provided, input files must be written
        to the test `workspace` prior to calling `run()`.
    check : callable, optional
        User defined function to evaluate results of the simulation.
        Takes `self` as input. This is a good place for assertions.
    plot: callable, optional
        User defined function to create plots of the results.
    compare: str or Comparison, optional
        Selects the comparison type or program name. If a string,
        must be a key into the `targets` dictionary, i.e. a valid
        program to use for comparison. Acceptable values: auto,
        mf6, mf6_regression, libmf6, mf2005, mfnwt, mflgr, mfnwt.
        If 'auto', the program to use is determined automatically
        by the contents of the comparison model/simulation folder.
    parallel : bool, optional
        Whether to test mf6 parallel capabilities.
    ncpus : int, optional
        Number of CPUs for mf6 parallel testing.
    htol : float, optional
        Tolerance for result comparisons.
    rclose : float, optional
        Residual tolerance for convergence
    overwrite : bool, optional
        Whether to overwrite existing output files in the workspace.
    verbose: bool, optional
        Whether to show verbose output
    xfail : bool, optional
        Whether the test is expected to fail

    """

    # tell pytest this class doesn't contain tests, don't collect it
    __test__ = False

    def __init__(
        self,
        name: str,
        workspace: str | os.PathLike,
        targets: dict[str, Path],
        api_func: Callable | None = None,
        build: Callable | None = None,
        check: Callable | None = None,
        plot: Callable | None = None,
        compare: str | Comparison | None = "auto",
        parallel: bool = False,
        ncpus: int = 1,
        htol: float | None = None,
        rclose: float | None = None,
        overwrite: bool = True,
        verbose: bool = False,
        xfail: bool | list[bool] = False,
        cargs: list | None = None,
    ):
        # make sure workspace exists
        workspace = Path(workspace).expanduser().absolute()
        assert workspace.is_dir(), f"{workspace} is not a valid directory"
        if verbose:
            print("Initializing test", name, "in workspace", workspace)

        self.name = name
        self.workspace = workspace
        self.targets = targets
        self.build = build
        self.check = check
        self.plot = plot
        self.parallel = parallel
        self.ncpus = [ncpus] if isinstance(ncpus, int) else ncpus
        self.api_func = api_func
        self.compare = Comparison(compare) if compare else None
        if self.compare in FP4_COMPARISONS:
            require_flopy4()
        self.outp = None
        self.htol = 0.001 if htol is None else htol
        self.rclose = 0.001 if rclose is None else rclose
        self.overwrite = overwrite
        self.verbose = verbose
        self.xfail = [xfail] if isinstance(xfail, bool) else xfail
        self.cargs = cargs

    def __repr__(self):
        return self.name

    # private

    def _compare_heads(
        self,
        cpth=None,
        extensions="hds",
        mf6=False,
        htol=0.001,
        cmp_dir="mf6_regression",
        workspace=None,
    ) -> bool:
        if isinstance(extensions, str):
            extensions = [extensions]

        if cpth:
            files1 = []
            files2 = []
            exfiles = []
            for file1 in self.outp:
                ext = os.path.splitext(file1)[1][1:]
                if ext.lower() in extensions:
                    # simulation file
                    pth = os.path.join(self.workspace, file1)
                    files1.append(pth)

                    # look for an exclusion file
                    pth = os.path.join(self.workspace, file1 + ".ex")
                    exfiles.append(pth if os.path.isfile(pth) else None)

                    # look for a comparison file
                    coutp = None
                    if mf6:
                        _, coutp = get_mf6_files(cpth / "mfsim.nam")
                    if coutp is not None:
                        for file2 in coutp:
                            ext = os.path.splitext(file2)[1][1:]
                            if ext.lower() in extensions:
                                files2.append(os.path.join(cpth, file2))
                    else:
                        files2.append(None)

            # todo: clean up namfile path detection?
            nf = next(iter(get_namefiles(cpth)), None)
            cmp_namefile = (
                None
                if self.compare
                in [Comparison.MF6, Comparison.MF6_REGRESSION, Comparison.LIBMF6]
                else os.path.basename(nf)
                if nf
                else None
            )
            if cmp_namefile is None:
                pth = None
            else:
                pth = os.path.join(cpth, cmp_namefile)

            for i in range(len(files1)):
                file1 = files1[i]
                ext = os.path.splitext(file1)[1][1:].lower()
                outfile = os.path.splitext(os.path.basename(file1))[0]
                outfile = os.path.join(self.workspace, outfile + "." + ext + ".cmp.out")
                file2 = None if files2 is None else files2[i]

                # set exfile
                exfile = None
                if file2 is None:
                    if len(exfiles) > 0:
                        exfile = exfiles[i]
                        if exfile is not None:
                            print(f"Exclusion file {i + 1}", os.path.basename(exfile))

                # make comparison
                success = compare_heads(
                    None,
                    pth,
                    precision="double",
                    text=EXTTEXT[ext],
                    outfile=outfile,
                    files1=file1,
                    files2=file2,
                    htol=htol,
                    difftol=True,
                    verbose=self.verbose,
                    exfile=exfile,
                )
                print(f"{EXTTEXT[ext]} comparison {i + 1}", self.name)
                if not success:
                    return False
            return True

        # otherwise it's a regression comparison
        ws = workspace or self.workspace
        files0, files1 = get_comparison_files(ws, extensions, cmp_dir)
        extension = "hds"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            outfile = os.path.splitext(os.path.basename(fpth0))[0]
            outfile = os.path.join(ws, outfile + f".{extension}.cmp.out")
            success = compare_heads(
                None,
                None,
                precision="double",
                htol=htol,
                text=EXTTEXT[extension],
                outfile=outfile,
                files1=fpth0,
                files2=fpth1,
                verbose=self.verbose,
            )
            print(
                f"{EXTTEXT[extension]} comparison {i + 1}"
                + f"{self.name} ({os.path.basename(fpth0)})"
            )
            if not success:
                return False
        return True

    def _compare_concentrations(
        self, extensions="ucn", htol=0.001, cmp_dir="mf6_regression", workspace=None
    ) -> bool:
        if isinstance(extensions, str):
            extensions = [extensions]

        ws = workspace or self.workspace
        files0, files1 = get_comparison_files(ws, extensions, cmp_dir)
        extension = "ucn"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            outfile = os.path.splitext(os.path.basename(fpth0))[0]
            outfile = os.path.join(ws, outfile + f".{extension}.cmp.out")
            success = compare_heads(
                None,
                None,
                precision="double",
                htol=htol,
                text=dependent_variable_text(fpth0, EXTTEXT[extension]),
                outfile=outfile,
                files1=fpth0,
                files2=fpth1,
                verbose=self.verbose,
            )
            print(
                (
                    f"{EXTTEXT[extension]} comparison {i + 1}"
                    + f"{self.name} ({os.path.basename(fpth0)})",
                )
            )
            if not success:
                return False
        return True

    def _compare_budgets(
        self, extensions="cbc", rclose=0.001, cmp_dir="mf6_regression", workspace=None
    ) -> bool:
        if isinstance(extensions, str):
            extensions = [extensions]
        ws = workspace or self.workspace
        files0, files1 = get_comparison_files(ws, extensions, cmp_dir)
        extension = "cbc"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            print(
                f"{EXTTEXT[extension]} comparison {i + 1}",
                f"{self.name} ({os.path.basename(fpth0)})",
            )
            # a budget file with nothing saved is written empty
            if os.path.getsize(fpth0) == 0 and os.path.getsize(fpth1) == 0:
                continue
            outname = os.path.splitext(os.path.basename(fpth0))[0]
            outfile = os.path.join(ws, f"{outname}.{extension}.cmp.out")
            success = compare_cell_budget(fpth0, fpth1, outfile=outfile, rclose=rclose)
            if not success:
                return False
        return True

    def _sim_dirs(self) -> list[Path]:
        """The test's MF6 simulation workspaces, relative to its own, in run order."""
        dirs = [
            Path(get_workspace(s)).absolute().relative_to(self.workspace)
            for s in self.sims
            if isinstance(s, MFSimulation)
        ]
        return dirs or [Path(".")]

    def _compare(self, comparison: Comparison):
        """
        Compare the main simulation's output with that of another simulation or model.

        compare : str
            The comparison executable name: mf6, mf6_regression, libmf6, mf2005,
            mfnwt, mflgr, or mfusg.
        """

        if self.verbose:
            print("Comparison test", self.name)

        # adjust htol if < IMS outer_dvclose, and rclose for budget comparisons
        htol = adjust_htol(self.workspace, self.htol)
        rclose = get_rclose(self.workspace)
        cmp_path = self.workspace / comparison.value
        if comparison == Comparison.MF6_REGRESSION:
            assert self._compare_heads(extensions=HDS_EXT, htol=htol), (
                "head comparison failed"
            )
            assert self._compare_budgets(extensions=CBC_EXT, rclose=rclose), (
                "budget comparison failed"
            )
            assert self._compare_concentrations(htol=htol), (
                "concentration comparison failed"
            )
        elif comparison in FP4_COMPARISONS:
            for sim_dir in self._sim_dirs():
                ws = self.workspace / sim_dir
                cmp_dir = os.path.relpath(
                    self.workspace / comparison.value / sim_dir, ws
                )
                htol = adjust_htol(ws, self.htol)
                rclose = get_rclose(ws)
                assert self._compare_heads(
                    extensions=HDS_EXT, htol=htol, cmp_dir=cmp_dir, workspace=ws
                ), f"head comparison failed: {sim_dir}"
                assert self._compare_budgets(
                    extensions=CBC_EXT, rclose=rclose, cmp_dir=cmp_dir, workspace=ws
                ), f"budget comparison failed: {sim_dir}"
                assert self._compare_concentrations(
                    htol=htol, cmp_dir=cmp_dir, workspace=ws
                ), f"concentration comparison failed: {sim_dir}"
        else:
            assert self._compare_heads(
                cpth=cmp_path,
                extensions=HDS_EXT,
                mf6=comparison in [Comparison.MF6, Comparison.LIBMF6],
                htol=htol,
            ), "head comparison failed"

    def _run(
        self,
        workspace: str | os.PathLike,
        target: str | os.PathLike,
        xfail: bool = False,
        cargs: str | None = None,
        ncpus: int = 1,
    ) -> tuple[bool, list[str]]:
        """
        Run a simulation or model with FloPy.

        workspace : str or path-like
            The simulation or model workspace
        target : str or path-like
            The target executable to use
        xfail : bool
            Whether to expect failure
        ncpus : int
            The number of CPUs for a parallel run
        """

        # make sure workspace exists
        workspace = Path(workspace).expanduser().absolute()
        assert workspace.is_dir(), f"Workspace not found: {workspace}"

        # make sure executable exists and framework knows about it.
        # shutil.which() on Windows (Python 3.12+) won't resolve a full
        # path whose extension isn't in PATHEXT (e.g. libmf6.dll), so
        # fall back to the target path itself when it doesn't resolve.
        resolved = shutil.which(str(target))
        tgt = Path(resolved) if resolved else Path(target)
        assert tgt.is_file(), f"Target executable not found: {target}"
        assert tgt in self.targets.values(), (
            "Targets must be explicitly registered with the test framework"
        )

        if self.verbose:
            print(f"Running {target} in {workspace}")

        # needed in _compare_heads()... todo: inject explicitly?
        nf = next(iter(get_namefiles(workspace)), None)
        self.cmp_namefile = (
            None
            if "mf6" in target.name or "libmf6" in target.name
            else os.path.basename(nf)
            if nf
            else None
        )

        # run the model
        try:
            # via MODFLOW API
            if "libmf6" in target.name and self.api_func:
                success, buff = self.api_func(target, workspace)
            # via MF6 executable
            elif "mf6" in target.name:
                # parallel test if configured
                if self.parallel and ncpus > 1:
                    print(f"Parallel test {self.name} on {self.ncpus} processes")
                    try:
                        success, buff = run_parallel(workspace, target, ncpus)
                    except Exception:
                        warn(
                            "MODFLOW 6 parallel test",
                            self.name,
                            f"failed with error:\n{format_exc()}",
                        )
                        success = False
                else:
                    # otherwise serial run
                    try:
                        success, buff = flopy.run_model(
                            target,
                            workspace / "mfsim.nam",
                            model_ws=workspace,
                            report=True,
                            cargs=cargs,
                        )
                    except Exception:
                        warn(
                            "MODFLOW 6 serial test",
                            self.name,
                            f"failed with error:\n{format_exc()}",
                        )
                        success = False
            else:
                # non-MF6 model
                try:
                    nf_ext = ".mpsim" if "mp7" in target.name else ".nam"
                    namefile = next(iter(workspace.glob(f"*{nf_ext}")), None)
                    assert namefile, f"Control file with extension {nf_ext} not found"
                    success, buff = flopy.run_model(
                        target, namefile, workspace, report=True
                    )
                except Exception:
                    warn(f"{target} model failed:\n{format_exc()}")
                    success = False

            if xfail:
                if success:
                    warn("MODFLOW 6 model should have failed!")
                    success = False
                else:
                    success = True

        except Exception:
            success = False
            warn(f"Unhandled error in comparison model {self.name}:\n{format_exc()}")

        return success, buff

    # public

    def run(self):
        """
        Run the test case.
        """

        # if build fn provided, build models/simulations and write input files
        if self.build:
            sims = self.build(self)
            sims = sims if isinstance(sims, Iterable) else [sims]
            sims = [sim for sim in sims if sim]  # filter Nones
            self.sims = sims
            nsims = len(sims)
            self.buffs = list(repeat(None, nsims))

            assert len(self.xfail) in [
                1,
                nsims,
            ], "Invalid xfail: expected a single boolean or one for each model"
            if len(self.xfail) == 1 and nsims:
                self.xfail = list(repeat(self.xfail[0], nsims))

            assert len(self.ncpus) in [
                1,
                nsims,
            ], "Invalid ncpus: expected a single integer or one for each model"
            if len(self.ncpus) == 1 and nsims:
                self.ncpus = list(repeat(self.ncpus[0], nsims))

            write_input(*sims, overwrite=self.overwrite, verbose=self.verbose)
        else:
            self.sims = [MFSimulation.load(sim_ws=self.workspace)]
            self.buffs = [None]
            assert len(self.xfail) == 1, "Invalid xfail: expected a single boolean"
            assert len(self.ncpus) == 1, "Invalid ncpus: expected a single integer"

        # run models/simulations
        for i, sim_or_model in enumerate(self.sims):
            tgts = self.targets
            workspace = get_workspace(sim_or_model)
            exe_path = (
                Path(sim_or_model.exe_name) if sim_or_model.exe_name else tgts["mf6"]
            )
            target = (
                exe_path
                if exe_path in tgts.values()
                else tgts.get(exe_path.stem, tgts["mf6"])
            )
            xfail = self.xfail[i]
            cargs = self.cargs[i] if self.cargs is not None else None
            ncpus = self.ncpus[i]
            success, buff = self._run(workspace, target, xfail, cargs, ncpus)
            self.buffs[i] = buff  # store model output for assertions later
            assert success, (
                f"{'Simulation' if 'mf6' in str(target) else 'Model'} "
                f"{'should have failed' if xfail else 'failed'}: {workspace}"
            )

        # setup and run comparison model(s), if enabled
        if self.compare:
            # get expected output files from main simulation
            _, self.outp = get_mf6_files(self.workspace / "mfsim.nam", self.verbose)

            # try to autodetect comparison type if enabled
            if self.compare == Comparison.AUTO:
                if self.verbose:
                    print("Auto-detecting comparison type")
                self.compare = detect_comparison(self.workspace)
            if self.compare:
                if self.verbose:
                    print(f"Using comparison type: {self.compare}")

                # copy simulation files to comparison workspace if mf6 regression
                if self.compare == Comparison.MF6_REGRESSION:
                    cmp_path = self.workspace / self.compare.value
                    if os.path.isdir(cmp_path):
                        if self.verbose:
                            print(f"Cleaning {cmp_path}")
                        shutil.rmtree(cmp_path)
                    if self.verbose:
                        print(
                            "Copying simulation files "
                            f"from {self.workspace} to {cmp_path}"
                        )
                    shutil.copytree(self.workspace, cmp_path)

                # run comparison simulation
                if self.compare in FP4_COMPARISONS:
                    # the simulation re-written by flopy4 with NetCDF input
                    workspace = self.workspace / self.compare.value
                    sim_dirs = self._sim_dirs()
                    write_fp4_netcdf(
                        self.workspace,
                        workspace,
                        self.compare.value.removeprefix("fp4_"),
                        sim_dirs,
                    )
                    for sim_dir in sim_dirs:
                        sim_ws = workspace / sim_dir
                        success, _ = self._run(sim_ws, self.targets["mf6"])
                        assert success, f"flopy4 NetCDF simulation failed: {sim_ws}"
                        check_fp4_netcdf_input(sim_ws)
                elif self.compare.value not in self.targets:
                    warn(
                        f"Couldn't find comparison program '{self.compare}', "
                        "skipping comparison"
                    )
                else:
                    # todo: don't hardcode workspace or assume agreement with
                    # test case simulation workspace, set & access simulation
                    # workspaces directly
                    workspace = self.workspace / self.compare.value
                    success, _ = self._run(
                        workspace,
                        self.targets.get(self.compare.value, self.targets["mf6"]),
                    )
                    assert success, f"Comparison model failed: {workspace}"

                # compare model results, if enabled
                if self.verbose and self.compare.value in self.targets:
                    print("Comparing results")
                self._compare(self.compare)

        # check results, if enabled
        if self.check:
            if self.verbose:
                print("Checking results")
            self.check(self)

        # plot results, if enabled
        if self.plot:
            if self.verbose:
                print("Plotting results")
            self.plot(self)

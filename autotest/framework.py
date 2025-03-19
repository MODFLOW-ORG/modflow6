import os
import shutil
from collections.abc import Iterable
from itertools import repeat
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen
from traceback import format_exc
from typing import Callable, Optional, Union
from warnings import warn

import flopy
import numpy as np
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
from flopy.utils.compare import compare_heads
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
        workspace: Union[str, os.PathLike],
        targets: dict[str, Path],
        api_func: Optional[Callable] = None,
        build: Optional[Callable] = None,
        check: Optional[Callable] = None,
        plot: Optional[Callable] = None,
        compare: Optional[str] = "auto",
        parallel=False,
        ncpus=1,
        htol=None,
        rclose=None,
        overwrite=True,
        verbose=False,
        xfail=False,
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
        self.outp = None
        self.htol = 0.001 if htol is None else htol
        self.rclose = 0.001 if rclose is None else rclose
        self.overwrite = overwrite
        self.verbose = verbose
        self.xfail = [xfail] if isinstance(xfail, bool) else xfail

    def __repr__(self):
        return self.name

    # private

    def _compare_heads(
        self, cpth=None, extensions="hds", mf6=False, htol=0.001
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
        files0, files1 = get_comparison_files(self.workspace, extensions)
        extension = "hds"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            outfile = os.path.splitext(os.path.basename(fpth0))[0]
            outfile = os.path.join(self.workspace, outfile + f".{extension}.cmp.out")
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

    def _compare_concentrations(self, extensions="ucn", htol=0.001) -> bool:
        if isinstance(extensions, str):
            extensions = [extensions]

        files0, files1 = get_comparison_files(self.workspace, extensions)
        extension = "ucn"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            outfile = os.path.splitext(os.path.basename(fpth0))[0]
            outfile = os.path.join(self.workspace, outfile + f".{extension}.cmp.out")
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
                (
                    f"{EXTTEXT[extension]} comparison {i + 1}"
                    + f"{self.name} ({os.path.basename(fpth0)})",
                )
            )
            if not success:
                return False
        return True

    def _compare_budgets(self, extensions="cbc", rclose=0.001) -> bool:
        if isinstance(extensions, str):
            extensions = [extensions]
        files0, files1 = get_comparison_files(self.workspace, extensions)
        extension = "cbc"
        for i, (fpth0, fpth1) in enumerate(zip(files0, files1)):
            print(
                f"{EXTTEXT[extension]} comparison {i + 1}",
                f"{self.name} ({os.path.basename(fpth0)})",
            )
            success = self._compare_budget_files(extension, fpth0, fpth1, rclose)
            if not success:
                return False
        return True

    def _compare_budget_files(self, extension, fpth0, fpth1, rclose=0.001) -> bool:
        success = True
        if os.stat(fpth0).st_size * os.stat(fpth0).st_size == 0:
            return success, ""
        outfile = os.path.splitext(os.path.basename(fpth0))[0]
        outfile = os.path.join(self.workspace, outfile + f".{extension}.cmp.out")
        fcmp = open(outfile, "w")
        fcmp.write("Performing CELL-BY-CELL to CELL-BY-CELL comparison\n")
        fcmp.write(f"{fpth0}\n")
        fcmp.write(f"{fpth1}\n\n")

        # open the files
        cbc0 = flopy.utils.CellBudgetFile(
            fpth0, precision="double", verbose=self.verbose
        )
        cbc1 = flopy.utils.CellBudgetFile(
            fpth1, precision="double", verbose=self.verbose
        )

        # build list of cbc data to retrieve
        avail0 = cbc0.get_unique_record_names()
        avail1 = cbc1.get_unique_record_names()
        avail0 = [t.decode().strip() for t in avail0]
        avail1 = [t.decode().strip() for t in avail1]

        # initialize list for storing totals for each budget term terms
        cbc_keys0 = []
        cbc_keys1 = []
        for t in avail0:
            t1 = t
            if t not in avail1:
                # check if RCHA or EVTA is available and use that instead
                # should be able to remove this once v6.3.0 is released
                if t[:-1] in avail1:
                    t1 = t[:-1]
                else:
                    raise Exception(f"Could not find {t} in {fpth1}")
            cbc_keys0.append(t)
            cbc_keys1.append(t1)

        # get list of times and kstpkper
        kk = cbc0.get_kstpkper()
        times = cbc0.get_times()

        # process data
        for key, key1 in zip(cbc_keys0, cbc_keys1):
            for idx, (k, t) in enumerate(zip(kk, times)):
                v0 = cbc0.get_data(kstpkper=k, text=key)[0]
                v1 = cbc1.get_data(kstpkper=k, text=key1)[0]
                if v0.dtype.names is not None:
                    v0 = v0["q"]
                    v1 = v1["q"]
                # skip empty vectors
                if v0.size < 1:
                    continue
                vmin = rclose
                if vmin < 1e-6:
                    vmin = 1e-6
                vmin_tol = 5.0 * vmin
                if v0.shape != v1.shape:
                    v0 = v0.flatten()
                    v1 = v1.flatten()
                idx = (abs(v0) > vmin) & (abs(v1) > vmin)
                diff = np.zeros(v0.shape, dtype=v0.dtype)
                diff[idx] = abs(v0[idx] - v1[idx])
                diffmax = diff.max()
                indices = np.where(diff == diffmax)[0]
                if diffmax > vmin_tol:
                    success = False
                    msg = (
                        f"{os.path.basename(fpth0)} - "
                        + f"{key:16s} "
                        + f"difference ({diffmax:10.4g}) "
                        + f"> {vmin_tol:10.4g} "
                        + f"at {indices.size} nodes "
                        + f" [first location ({indices[0] + 1})] "
                        + f"at time {t} "
                    )
                    fcmp.write(f"{msg}\n")
                    if self.verbose:
                        print(msg)

        fcmp.close()
        return success

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
        else:
            assert self._compare_heads(
                cpth=cmp_path,
                extensions=HDS_EXT,
                mf6=comparison in [Comparison.MF6, Comparison.LIBMF6],
                htol=htol,
            ), "head comparison failed"

    def _run(
        self,
        workspace: Union[str, os.PathLike],
        target: Union[str, os.PathLike],
        xfail: bool = False,
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

        # make sure executable exists and framework knows about it
        tgt = Path(shutil.which(target))
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
        Run the test case end-to-end.

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
            ncpus = self.ncpus[i]
            success, buff = self._run(workspace, target, xfail, ncpus)
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

                # run comparison simulation if libmf6 or mf6 regression
                if self.compare in [Comparison.MF6_REGRESSION, Comparison.LIBMF6]:
                    if self.compare.value not in self.targets:
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

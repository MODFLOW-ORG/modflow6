"""
Release times after the end of the simulation are skipped
unless EXTEND_TRACKING is on. Check through the API that a
skipped release leaves no unfilled slots in the committed
particle store. The release scheduler hands such times to
the last time step, so the store is checked after each step.

Two cases: one without EXTEND_TRACKING (late release time is
skipped) and one with it (late release time is used).
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtskip"
cases = [simname, f"{simname}ext"]
nstp = 3
releasetimes = [0.0, 5.0]  # second is after the simulation ends


def build_models(idx, test):
    name = test.name
    sim = FlopyReadmeCase.get_gwf_sim(name, test.workspace, test.targets["libmf6"])
    sim.tdis.perioddata = [(FlopyReadmeCase.perlen, nstp, FlopyReadmeCase.tsmult)]

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdis(
        prt,
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )
    flopy.mf6.ModflowPrtmip(prt, porosity=FlopyReadmeCase.porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        nreleasepts=len(FlopyReadmeCase.releasepts_prt),
        packagedata=FlopyReadmeCase.releasepts_prt,
        nreleasetimes=len(releasetimes),
        releasetimes=[[t] for t in releasetimes],
        extend_tracking="ext" in name,
    )
    flopy.mf6.ModflowPrtoc(prt, trackcsv_filerecord=[f"{prt_name}.trk.csv"])
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(name, "gwf"),
        exgmnameb=prt_name,
    )
    ems = flopy.mf6.ModflowEms(sim, filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx]
    model_ws = Path(model_ws or ".")
    prt_name = get_model_name(name, "prt").upper()
    nrpt = len(FlopyReadmeCase.releasepts_prt)
    extend = "ext" in name

    mf6 = ModflowApi(exe, working_directory=model_ws)
    mf6.initialize()
    try:
        kstp = 0
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.update()
            kstp += 1

            def get(var):
                return mf6.get_value(mf6.get_var_address(var, prt_name, "PRP1"))

            irpt = get("PLIRPT")
            trelease = get("PLTRELEASE")

            # the late release time is only used with extend tracking,
            # and only on the last step, which the scheduler gives it to
            ntimes = 2 if extend and kstp == nstp else 1
            assert irpt.size == nrpt * ntimes, (
                f"step {kstp}: expected {nrpt * ntimes} particles, got {irpt.size}"
            )
            assert np.all((irpt >= 1) & (irpt <= nrpt))
            assert np.allclose(np.unique(trelease), releasetimes[:ntimes])
            assert get("NPARTICLES")[0] == irpt.size
    finally:
        mf6.finalize()

    return True, open(model_ws / "mfsim.stdout").readlines()


def check_output(idx, test):
    stdout = (test.workspace / "mfsim.stdout").read_text()
    nwarn = stdout.count("Skipping release time falling after the end")
    assert nwarn == (0 if "ext" in test.name else 1)


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()

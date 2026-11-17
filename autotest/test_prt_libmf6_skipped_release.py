"""
Release times that won't be used are dropped when the PRP
package is loaded: negative times, and times after the end of
the simulation unless EXTEND_TRACKING is on. Check through the
API that the committed particle store holds only the particles
actually released, after each time step. The release scheduler
hands times after the end to the last time step, so a dropped
time would otherwise leave unfilled slots there.

Cases:
- late:     one time after the end, dropped
- lateext:  one time after the end, kept (EXTEND_TRACKING)
- alllate:  all times after the end, dropped; this must not
            trigger the default release at t=0
- neg:      one negative time, dropped
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from prt_test_utils import FlopyReadmeCase, get_model_name

nstp = 3
perlen = FlopyReadmeCase.perlen  # simulation ends at t=1

LATE_WARNING = "Skipping release time falling after the end"
NEG_WARNING = "Skipping negative release time"

# releasetimes: explicit release times in the PRP input
# extend: whether EXTEND_TRACKING is on
# released: release times expected in the store after each step
# warnings: expected count of each warning in mfsim.stdout
cases = {
    "prtlate": {
        "releasetimes": [0.0, 5.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {LATE_WARNING: 1, NEG_WARNING: 0},
    },
    "prtlateext": {
        "releasetimes": [0.0, 5.0],
        "extend": True,
        "released": [[0.0], [0.0], [0.0, 5.0]],
        "warnings": {LATE_WARNING: 0, NEG_WARNING: 0},
    },
    "prtalllate": {
        "releasetimes": [5.0, 6.0],
        "extend": False,
        "released": [[], [], []],
        "warnings": {LATE_WARNING: 2, NEG_WARNING: 0},
    },
    "prtneg": {
        "releasetimes": [-1.0, 0.0],
        "extend": False,
        "released": [[0.0], [0.0], [0.0]],
        "warnings": {LATE_WARNING: 0, NEG_WARNING: 1},
    },
}
names = list(cases)


def build_models(idx, test):
    name = test.name
    case = cases[name]
    sim = FlopyReadmeCase.get_gwf_sim(name, test.workspace, test.targets["libmf6"])
    sim.tdis.perioddata = [(perlen, nstp, FlopyReadmeCase.tsmult)]

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
        nreleasetimes=len(case["releasetimes"]),
        releasetimes=[[t] for t in case["releasetimes"]],
        extend_tracking=case["extend"],
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

    name = names[idx]
    case = cases[name]
    model_ws = Path(model_ws or ".")
    prt_name = get_model_name(name, "prt").upper()
    nrpt = len(FlopyReadmeCase.releasepts_prt)

    mf6 = ModflowApi(exe, working_directory=model_ws)
    mf6.initialize()
    try:
        kstp = 0
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.update()
            kstp += 1

            def get(var):
                return mf6.get_value(mf6.get_var_address(var, prt_name, "PRP1"))

            released = case["released"][kstp - 1]
            irpt = get("PLIRPT")
            trelease = get("PLTRELEASE")
            assert irpt.size == nrpt * len(released), (
                f"step {kstp}: expected {nrpt * len(released)} particles, "
                f"got {irpt.size}"
            )
            assert np.all((irpt >= 1) & (irpt <= nrpt))
            assert np.allclose(np.unique(trelease), released)
            assert get("NPARTICLES")[0] == irpt.size
        assert kstp == nstp
    finally:
        mf6.finalize()

    return True, open(model_ws / "mfsim.stdout").readlines()


def check_output(idx, test):
    stdout = (test.workspace / "mfsim.stdout").read_text()
    for warning, count in cases[test.name]["warnings"].items():
        assert stdout.count(warning) == count, warning


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(names))
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

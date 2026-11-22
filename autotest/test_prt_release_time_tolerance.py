"""
Test the PRP release time coincidence tolerance.

The default RELEASE_TIME_TOLERANCE scales with the simulation
length T (the sum of all stress period lengths), so whether two
nearby release times are merged depends on their spacing relative
to T, not on the model's time units. Separately, release times
indistinguishable at floating-point precision are always merged,
even if the specified tolerance is below the ulp spacing.

Each case releases particles from one point at two explicit times
(or one explicit time and one period-block time), and checks how
many distinct releases occurred.

Cases:
    - floor: T = 2e10, period-block FIRST release in period 2 at
      t = 1e10, and an explicit release 10 ulps earlier. Tolerance
      is set far below the ulp spacing. Expect a single release.
    - sclsml: T = 1, explicit releases 1e-6 * T apart. Expect two.
      (The old fixed default, ~2.2e-5, merged these.)
    - sclbig: T = 1e6, explicit releases 1e-6 * T apart. Expect two.
    - short: T = 1, explicit releases 1e-7 apart. Expect two.
      (The old fixed default merged these.)
    - long: T = 1e6, explicit releases 1e-4 apart. Expect one.
      (The old fixed default kept these separate.)
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prttol"

# case name: (period lengths, explicit release times, period block, rttol, nrel)
FLOOR_T = 1e10
EPS = np.finfo(float).eps
DEFAULT_FRAC = 100 * EPS * 1e5
CASES = {
    "floor": (
        [FLOOR_T, FLOOR_T],
        [FLOOR_T - 10 * np.spacing(FLOOR_T)],
        {0: [], 1: [("FIRST",)]},
        1e-12,
        1,
    ),
    "sclsml": ([1.0], [0.5, 0.5 + 1e-6], None, None, 2),
    "sclbig": ([1e6], [0.5e6, 0.5e6 + 1.0], None, None, 2),
    "short": ([1.0], [0.5, 0.5 + 1e-7], None, None, 2),
    "long": ([1e6], [0.5e6, 0.5e6 + 1e-4], None, None, 1),
}
cases = [f"{simname}{k}" for k in CASES]


def get_case(name):
    return CASES[name.removeprefix(simname)]


def set_tdis(sim, perlens):
    tdis = sim.get_package("tdis")
    tdis.nper = len(perlens)
    tdis.perioddata = [(p, 1, 1.0) for p in perlens]


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    perlens, releasetimes, perioddata, rttol, _ = get_case(name)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, version="mf6", sim_ws=prt_ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=len(perlens),
        perioddata=[(p, 1, 1.0) for p in perlens],
    )
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # write release times with full precision
    prt_ws.mkdir(parents=True, exist_ok=True)
    releasetimes_path = prt_ws / "releasetimes.txt"
    with open(releasetimes_path, "w") as f:
        for t in releasetimes:
            f.write(f"{float(t)!r}\n")

    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=1,
        packagedata=FlopyReadmeCase.releasepts_prt[:1],
        nreleasetimes=len(releasetimes),
        releasetimes=f"open/close {releasetimes_path.name}",
        perioddata=perioddata,
        release_time_tolerance=rttol,
        trackcsv_filerecord=[f"{prt_name}.prp.trk.csv"],
    )
    flopy.mf6.ModflowPrtoc(prt, pname="oc", trackcsv_filerecord=[f"{prt_name}.trk.csv"])
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_ws / f"{gwf_name}.hds"),
            ("GWFBUDGET", gwf_ws / f"{gwf_name}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(test):
    perlens = get_case(test.name)[0]
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )
    set_tdis(gwf_sim, perlens)
    prt_sim = build_prt_sim(
        test.name, test.workspace, test.workspace / "prt", test.targets["mf6"]
    )
    return gwf_sim, prt_sim


def check_output(test):
    name = test.name
    perlens, _, _, rttol, nrel = get_case(name)
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(name, "prt")

    # check the number of distinct releases
    pls = pd.read_csv(prt_ws / f"{prt_name}.trk.csv")
    releases = pls[pls["ireason"] == 0]
    assert len(releases) == nrel, releases
    assert releases["trelease"].nunique() == nrel

    # check the effective tolerance is reported in the listing file
    lines = [l.strip() for l in open(prt_ws / f"{prt_name}.lst")]
    tol_lines = [l for l in lines if l.startswith("RELEASE TIME TOLERANCE =")]
    assert len(tol_lines) == 1
    tol_value = float(tol_lines[0].split("=")[1].split("(")[0])
    if rttol is None:
        assert "(DEFAULT)" in tol_lines[0]
        assert np.isclose(tol_value, DEFAULT_FRAC * sum(perlens), rtol=1e-6)
    else:
        assert "(SPECIFIED)" in tol_lines[0]
        assert np.isclose(tol_value, rttol, rtol=1e-6)


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t),
        check=lambda t: check_output(t),
        targets=targets,
        compare=None,
    )
    test.run()

"""
Regression test for the transient (STORAGE) SFR reach-aquifer exchange, issue
#2884.

With the SFR STORAGE option (transient kinematic-wave channel routing), a dry
gaining reach (aquifer head above the streambed top) with no inflow used to stay
dry with exactly zero groundwater discharge. At a small time step, zero reach
depth is a stable fixed point of the reach solver because the streambed
conductance and saturation vanish at zero depth, so the reach never rewets from
groundwater (it does rewet with the STORAGE option off, the steady solver). The
transient depth of a gaining reach is now found by a bracketed bisection, so it
rewets for any time-step size.

The model is a five-reach rectangular network with the aquifer held one meter
above the streambed by CHD and no inflow anywhere. It is run with STORAGE on at a
time step near a Courant number of one (the regime the bug was worst in). The
test checks that every reach rewets (stage above the streambed top) and that the
network receives groundwater, rather than staying trapped at zero discharge.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ("sfr-drygain",)

NREACH = 5
DELR = 100.0
SLOPE = 0.001
WIDTH = 5.0
HK = 0.1
GAIN_ABOVE = 1.0
RTP = [10.0 - SLOPE * ((j + 0.5) * DELR) for j in range(NREACH)]
# time step near Courant ~ 1 (celerity ~1 m/s over a 100 m reach) with enough
# steps to reach quasi-steady state (the exchange has fully settled well before
# this; the unfixed solver stays trapped at zero depth for any step count)
DELT, NSTP = 100.0, 25


def build_models(idx, test):
    name = cases[idx]
    head = [RTP[j] + GAIN_ABOVE for j in range(NREACH)]

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim, time_units="seconds", nper=1, perioddata=[(DELT * NSTP, NSTP, 1.0)]
    )
    flopy.mf6.ModflowIms(
        sim,
        complexity="complex",
        outer_maximum=500,
        inner_maximum=200,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=2,
        nrow=1,
        ncol=NREACH,
        delr=DELR,
        delc=100.0,
        top=12.0,
        botm=[4.0, 0.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=np.array([head, head])[:, None, :])
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=[1, 0], k=10.0, k33=1000.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=[1, 0], ss=1e-5, sy=0.2, transient={0: True})
    # hold the aquifer head above the streambed with CHD in the lower layer
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(1, 0, j), float(head[j])] for j in range(NREACH)],
    )

    packagedata, connectiondata = [], []
    for j in range(NREACH):
        ncon = 1 if j in (0, NREACH - 1) else 2
        packagedata.append(
            (j, (0, 0, j), DELR, WIDTH, SLOPE, RTP[j], 1.0, HK, 0.03, ncon, 1.0, 0)
        )
        conn = [j]
        if j > 0:
            conn.append(j - 1)
        if j < NREACH - 1:
            conn.append(-(j + 1))
        connectiondata.append(conn)

    sfr_obs = {
        f"{name}.sfr.obs.csv": [(f"stage{j + 1}", "stage", (j,)) for j in range(NREACH)]
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        storage=True,
        nreaches=NREACH,
        packagedata=packagedata,
        connectiondata=connectiondata,
        budget_filerecord=f"{name}.sfr.cbc",
        perioddata=[],  # no inflow anywhere
    )
    sfr.obs.initialize(filename=f"{name}.sfr.obs", continuous=sfr_obs, print_input=True)
    flopy.mf6.ModflowGwfoc(
        gwf, budget_filerecord=f"{name}.cbc", saverecord=[("BUDGET", "ALL")]
    )
    return sim


def check_output(idx, test):
    name = cases[idx]
    strtop = np.array(RTP)

    # reach stage at the final (quasi-steady) time step
    obs = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    stage = np.array([obs[f"STAGE{j + 1}"][-1] for j in range(NREACH)])
    depth = stage - strtop
    print("reach depth:", depth)

    # every gaining reach rewets instead of staying trapped at zero depth
    assert np.all(depth > 1e-6), (
        f"a gaining reach did not rewet (depth {depth}); the reach is trapped "
        "at zero groundwater discharge"
    )

    # the network receives groundwater (net exchange is materially nonzero)
    cbc = flopy.utils.CellBudgetFile(
        test.workspace / f"{name}.sfr.cbc", precision="double"
    )
    gwf = cbc.get_data(text="GWF", kstpkper=cbc.get_kstpkper()[-1])[0]
    qgwf = np.array([r["q"] for r in gwf])
    print("reach qgwf:", qgwf)
    assert abs(float(qgwf.sum())) > 1e-3, (
        f"network groundwater exchange is ~zero (sum {qgwf.sum()}); the reaches "
        "are trapped"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()

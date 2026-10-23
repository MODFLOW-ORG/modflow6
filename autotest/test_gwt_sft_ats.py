"""
Tests that SFT gives the same result when ATS repeats a time step as when the
shorter time step is taken from the start.

Three reaches without SFR reach storage carry an inflow that rises from 10 to
100 m3/d, and from 50 to 500 g/m3, at the start of the second stress period;
one outer iteration fails the first, 1 d step of that period, and ATS repeats
it with a fifth of its length.

Cases:
  - sft_ats : the reach concentrations match those of a simulation without ATS
              that takes the same time steps.
"""

import flopy
import numpy as np
import pytest
from compare import Comparison
from framework import TestFramework

nreach = 3
rlen, rwid, rgrd, man = 1000.0, 10.0, 1.0e-4, 0.1
c0 = 50.0
c_in = 500.0
# the first period and the time steps that ATS takes in the second
perlen1 = 1.0
perlen2 = 2.0
dt_retry = 0.2


def build_sim(ws, exe, periods, ats=None):
    """Stream reaches over an aquifer that they do not exchange water with.

    Each period is (length, number of time steps, inflow, concentration of the
    inflow); ats is the ATS period data.
    """
    sim = flopy.mf6.MFSimulation(sim_name="sftats", sim_ws=ws, exe_name=exe)
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(p[0], p[1], 1.0) for p in periods],
    )
    if ats is not None:
        tdis.ats.initialize(maxats=len(ats), perioddata=ats, filename="sftats.ats")

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    sim.register_ims_package(
        flopy.mf6.ModflowIms(
            sim,
            outer_dvclose=1e-10,
            inner_dvclose=1e-12,
            linear_acceleration="bicgstab",
            filename="gwf.ims",
        ),
        ["gwf"],
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=nreach, delr=rlen, delc=100.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1.0e-6, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, j), 5.0] for j in range(nreach)]
    )
    # the streambed is impermeable, so the reaches exchange no water with the
    # aquifer, and the SFR Package does not route water through reach storage
    flopy.mf6.ModflowGwfsfr(
        gwf,
        length_conversion=1.0,
        time_conversion=86400.0,
        nreaches=nreach,
        packagedata=[
            (
                i,
                (0, 0, i),
                rlen,
                rwid,
                rgrd,
                9.0,
                1.0,
                0.0,
                man,
                1 if i in (0, 2) else 2,
                1.0,
                0,
            )
            for i in range(nreach)
        ],
        connectiondata=[[0, -1], [1, 0, -2], [2, 1]],
        perioddata={k: [(0, "inflow", p[2])] for k, p in enumerate(periods)},
        pname="SFR-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord="gwf.cbc",
        saverecord=[("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    # one outer iteration converges only if the concentration changes by less
    # than 60 g/m3
    sim.register_ims_package(
        flopy.mf6.ModflowIms(
            sim,
            outer_maximum=1,
            outer_dvclose=60.0,
            inner_dvclose=1e-12,
            linear_acceleration="bicgstab",
            filename="gwt.ims",
        ),
        ["gwt"],
    )
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=nreach, delr=rlen, delc=100.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtsft(
        gwt,
        flow_package_name="SFR-1",
        packagedata=[(i, c0) for i in range(nreach)],
        reachperioddata={k: [(0, "inflow", p[3])] for k, p in enumerate(periods)},
        concentration_filerecord="gwt.sft.ucn",
        budget_filerecord="gwt.sft.cbc",
        pname="SFT-1",
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def reach_conc(ws):
    f = flopy.utils.HeadFile(
        ws / "gwt.sft.ucn", text="CONCENTRATION", precision="double"
    )
    times = np.array(f.get_times())
    return times, np.array([f.get_data(totim=t).ravel() for t in times])


def build_models(idx, test):
    # the first period is at steady state; ATS starts the second with a 1 d
    # time step and repeats it with a fifth of that length
    periods = [(perlen1, 1, 10.0, c0), (perlen2, 1, 100.0, c_in)]
    ats = [(1, 1.0, dt_retry, 1.0, 1.0, 5.0)]
    return build_sim(test.workspace, "mf6", periods, ats=ats)


def check_output(idx, test):
    ws = test.workspace
    listing = (ws / "mfsim.lst").read_text()
    assert "Failed solution for step" in listing, "no time step was repeated"
    times, conc = reach_conc(ws)

    # the same time steps without ATS, as stress periods of one step each
    steps = np.diff(times)
    assert np.allclose(steps[1:], dt_retry), f"time steps {np.diff(times)}"
    periods = [(perlen1, 1, 10.0, c0)] + [(dt, 1, 100.0, c_in) for dt in steps]
    ref_ws = ws / "reference"
    sim = build_sim(ref_ws, str(test.targets["mf6"]), periods)
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"reference simulation failed\n{buff}"
    ref_times, ref_conc = reach_conc(ref_ws)
    assert np.allclose(ref_times, times), f"times {times} and {ref_times}"
    assert np.allclose(conc, ref_conc, rtol=1e-9, atol=1e-9), (
        f"reach concentrations with ATS\n{conc}\nwithout\n{ref_conc}"
    )


@pytest.mark.parametrize(
    "compare",
    [
        "auto",
        pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4),
    ],
)
def test_mf6model(compare, function_tmpdir, targets):
    test = TestFramework(
        name="sft_ats",
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(0, t),
        check=lambda t: check_output(0, t),
        compare=compare,
    )
    test.run()

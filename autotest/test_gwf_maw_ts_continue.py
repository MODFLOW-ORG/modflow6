"""A PERIOD RATE, WELL_HEAD, or AUXILIARY value linked to a time series
keeps tracking that time series in later periods whose own PERIOD block
doesn't reissue it, for GWF-MAW.

Model geometry: 1-layer, 1-row, 3-column; CHD at flanking columns; single
    MAW well in the centre column; 3 steady-state stress periods, 1 step
    each.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "maw"
cases = ["maw_ratecont", "maw_wellheadcont", "maw_auxcont"]

nper = 3
perlen = 1.0
period_data = [(perlen, 1, 1.0)] * nper

nlay, nrow, ncol = 1, 1, 3
delr = delc = 100.0
top = 200.0
botm = [0.0]
k11 = 1.0
chd_head = 100.0

radius = 0.1
well_bot = -10.0
strt = 90.0

rate_vals = [-1.0e-3, -2.0e-3, -3.0e-3]
well_head_vals = [80.0, 60.0, 40.0]
aux_vals = [12.0, 24.0, 36.0]

# TS times at period-end; LINEAREND gives exact per-period values.
ts_times = [0.0, 1.0, 2.0, 3.0]
ts_rate = [rate_vals[0]] + rate_vals
ts_head = [well_head_vals[0]] + well_head_vals
ts_aux = [aux_vals[0]] + aux_vals


def _base_sim(ws, name):
    """Return a minimal GWF sim ready for a MAW package to be attached."""
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=period_data)
    flopy.mf6.ModflowIms(sim, print_option="NONE")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=k11)
    flopy.mf6.ModflowGwfic(gwf, strt=chd_head)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), chd_head], [(0, 0, 2), chd_head]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )
    return sim, gwf


def _get_rate_continue(ws, name):
    """RATE linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "active"), (0, "rate", "rate_ts")],
            1: [(0, "status", "active")],
            2: [(0, "status", "active")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_rate", "RATE", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_rate)),
        time_series_namerecord=["rate_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_wellhead_continue(ws, name):
    """WELL_HEAD linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "constant"), (0, "well_head", "head_ts")],
            1: [(0, "status", "constant")],
            2: [(0, "status", "constant")],
        },
        pname="maw-1",
    )
    maw.obs.initialize(
        filename=f"{name}.{paktest}.obs",
        continuous={f"{name}.{paktest}.obs.csv": [("well_head", "HEAD", (0,))]},
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_head)),
        time_series_namerecord=["head_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


def _get_aux_continue(ws, name):
    """AUXILIARY linked to a time series; periods 1/2 reissue STATUS only."""
    sim, gwf = _base_sim(ws, name)
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        auxiliary=["conc"],
        budget_filerecord=f"{name}.{paktest}.cbc",
        print_input=True,
        packagedata=[(0, radius, well_bot, strt, "THIEM", 1, 0.0)],
        connectiondata=[(0, 0, (0, 0, 1), top, well_bot, 0.0, 0.0)],
        perioddata={
            0: [(0, "status", "active"), (0, "AUXILIARY", "conc", "conc_ts")],
            1: [(0, "status", "active")],
            2: [(0, "status", "active")],
        },
        pname="maw-1",
    )
    maw.ts.initialize(
        filename=f"{name}.{paktest}.ts",
        timeseries=list(zip(ts_times, ts_aux)),
        time_series_namerecord=["conc_ts"],
        interpolation_methodrecord=["linearend"],
    )
    return sim


_builders = [_get_rate_continue, _get_wellhead_continue, _get_aux_continue]


def build_models(idx, test):
    name = cases[idx]
    return _builders[idx](test.workspace, name)


def check_output(idx, test):
    name = cases[idx]
    if idx == 0:
        fname = test.workspace / f"{name}.{paktest}.obs.csv"
        tc = np.genfromtxt(fname, names=True, delimiter=",")
        obs = tc["WELL_RATE"]
        expected = rate_vals
        label = "RATE"
    elif idx == 1:
        fname = test.workspace / f"{name}.{paktest}.obs.csv"
        tc = np.genfromtxt(fname, names=True, delimiter=",")
        obs = tc["WELL_HEAD"]
        expected = well_head_vals
        label = "WELL_HEAD"
    else:
        fname = test.workspace / f"{name}.{paktest}.cbc"
        cbc = flopy.utils.CellBudgetFile(fname, precision="double")
        data = cbc.get_data(text="AUXILIARY")
        auxcol = data[0].dtype.names[-1]
        obs = np.array([rec[auxcol][0] for rec in data])
        expected = aux_vals
        label = "AUXILIARY"

    assert len(obs) == nper, f"Expected {nper} {label} records, got {len(obs)}"
    assert np.allclose(obs, expected), (
        f"{label} expected to track its time series per period {expected} "
        f"(periods 1/2's PERIOD block doesn't reissue it), got {obs}."
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

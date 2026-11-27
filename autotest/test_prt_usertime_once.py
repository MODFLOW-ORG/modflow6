"""
User tracking times coinciding with other particle tracking events.

The same instant can be a particle's release time, a cell or subcell exit
time, a time step boundary and a user tracking time all at once. PRT should
report each event at the earliest opportunity, and each tracking time at
most once per particle. Events of different kinds at the same time are all
reported.

The grid is rectilinear, so it can be built as both DIS and DISV. Cases,
each on DIS with Pollock's method and on DISV with Pollock's method and
with the ternary method:

    face:     release on a cell face, with a tracking time equal to the
              release time.
    corner:   the same, released on a cell corner.
    exit:     tracking times equal to cell exit times, taken from a run
              without tracking times.
    tstep:    tracking times equal to time step ends.

Two cases use the STAY dry tracking method. In `stay`, particles sit in a
dry cell. In `arrival`, a particle moves into a dry cell partway through
the time step.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import get_model_name, get_track_dtype

simname = "ties"
nrow = ncol = 10
top, botm = 1.0, [0.0]
porosity = 0.1
perlen, nstp = 1.0, 4  # step ends at 0.25, 0.5 and 0.75 are exact

# release cell (row 2, column 2) spans x in [2, 3], y in [7, 8].
# flow is from the top left corner to the bottom right, i.e. +x, -y,
# so a particle on the east face or south east corner exits at once.
RELEASE_CELL = (0, 2, 2)
RELEASE_POINTS = {
    "face": (3.0, 7.5, 0.5),
    "corner": (3.0, 7.0, 0.5),
    # off the cell center and the grid's line of symmetry (x + y = 10),
    # along which the particle would track exactly on an internal edge
    # between the ternary method's triangular subcells
    "exit": (2.4, 7.3, 0.5),
    "tstep": (2.4, 7.3, 0.5),
}
GRIDS = {"dis": "pollock", "disv": "pollock", "disvt": "ternary"}
CODES = {"face": "f", "corner": "c", "exit": "x", "tstep": "s"}
cases = [f"{simname}{CODES[case]}_{grid}" for case in RELEASE_POINTS for grid in GRIDS]

RELEASE, EXIT, TIMESTEP, TERMINATE, USERTIME = 0, 1, 2, 3, 5


def parse(name):
    code, grid = name[len(simname) :].split("_")
    case = next(k for k, v in CODES.items() if v == code)
    return case, grid


def float_exactly(sim):
    # write floats with enough digits to round-trip, since tracking
    # times must equal computed exit times to the last bit
    sim.simulation_data.float_precision = 17
    sim.simulation_data.float_characters = 26


def add_dis(model, grid):
    if grid == "dis":
        cls = (
            flopy.mf6.ModflowGwfdis
            if isinstance(model, flopy.mf6.ModflowGwf)
            else flopy.mf6.ModflowPrtdis
        )
        cls(
            model,
            nlay=1,
            nrow=nrow,
            ncol=ncol,
            delr=1.0,
            delc=1.0,
            top=top,
            botm=botm,
        )
    else:
        cls = (
            flopy.mf6.ModflowGwfdisv
            if isinstance(model, flopy.mf6.ModflowGwf)
            else flopy.mf6.ModflowPrtdisv
        )
        cls(
            model,
            **get_disv_kwargs(
                nlay=1,
                nrow=nrow,
                ncol=ncol,
                delr=1.0,
                delc=1.0,
                tp=top,
                botm=botm,
            ),
        )


def cellid(grid, kij):
    k, i, j = kij
    return (k, i, j) if grid == "dis" else (k, i * ncol + j)


def build_gwf_sim(name, ws, mf6, grid):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=ws)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )
    gwfname = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    add_dis(gwf, grid)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [cellid(grid, (0, 0, 0)), top],
            [cellid(grid, (0, nrow - 1, ncol - 1)), 0.0],
        ],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(sim)
    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6, grid, case, tracktimes):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=prt_ws)
    float_exactly(sim)
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )
    prtname = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    add_dis(prt, grid)
    flopy.mf6.ModflowPrtmip(prt, porosity=porosity)
    x, y, z = RELEASE_POINTS[case]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=1,
        packagedata=[[0, cellid(grid, RELEASE_CELL), x, y, z]],
        perioddata={0: [("FIRST",)]},
        dev_forceternary=GRIDS[grid] == "ternary",
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[f"{prtname}.trk"],
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
        ntracktimes=len(tracktimes) if tracktimes else None,
        tracktimes=[(t,) for t in tracktimes] if tracktimes else None,
    )
    gwfname = get_model_name(name, "gwf")
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", Path(gwf_ws) / f"{gwfname}.hds"),
            ("GWFBUDGET", Path(gwf_ws) / f"{gwfname}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prtname}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def read_track(ws, name):
    prtname = get_model_name(name, "prt")
    ws = Path(ws)
    dt = get_track_dtype(ws / f"{prtname}.trk.hdr")
    return pd.DataFrame(np.fromfile(ws / f"{prtname}.trk", dtype=dt))


def exit_times(name, gwf_sim, prt_ws, mf6, grid):
    """Run without tracking times, and return the first few cell exit times."""
    gwf_sim.write_simulation(silent=True)
    ok, _ = gwf_sim.run_simulation(silent=True)
    assert ok
    pre_ws = Path(prt_ws).parent / "prt_pre"
    sim = build_prt_sim(name, gwf_sim.sim_path, pre_ws, mf6, grid, "exit", None)
    sim.write_simulation(silent=True)
    ok, _ = sim.run_simulation(silent=True)
    assert ok
    trk = read_track(pre_ws, name)
    times = trk[(trk.ireason == EXIT) & (trk.t > 0)].t.unique()
    assert len(times) >= 3
    return [float(t) for t in times[:3]]


def build_models(test, grid, case):
    mf6 = test.targets["mf6"]
    gwf_sim = build_gwf_sim(test.name, test.workspace, mf6, grid)
    prt_ws = test.workspace / "prt"
    if case in ("face", "corner"):
        tracktimes = [0.0, 0.1, 0.2]
    elif case == "exit":
        tracktimes = exit_times(test.name, gwf_sim, prt_ws, mf6, grid)
    elif case == "tstep":
        tracktimes = [k * perlen / nstp for k in range(1, nstp)]
    test.tracktimes = tracktimes
    prt_sim = build_prt_sim(
        test.name, test.workspace, prt_ws, mf6, grid, case, tracktimes
    )
    return gwf_sim, prt_sim


def check_usertimes_once(trk, tracktimes):
    """Each tracking time is reported once."""
    ut = trk[trk.ireason == USERTIME]
    assert sorted(ut.t.tolist()) == sorted(tracktimes), (
        f"expected each of {tracktimes} once, got {ut.t.tolist()}"
    )
    # time never goes backwards, and no exit time is negative
    assert (trk.t.diff().dropna() >= 0).all(), trk
    assert trk.iloc[0].ireason == RELEASE


def check_output(test, grid, case):
    trk = read_track(test.workspace / "prt", test.name)
    tracktimes = test.tracktimes
    check_usertimes_once(trk, tracktimes)

    if case in ("face", "corner"):
        # release, then the tracking time at the release time, then exit
        assert trk.iloc[1].ireason == USERTIME and trk.iloc[1].t == 0.0
        assert trk.iloc[2].ireason == EXIT
        assert trk.iloc[1].icell == trk.iloc[2].icell == trk.iloc[0].icell
        # on the face it exits through, so (all but) immediately. ternary
        # nudges barycentric coordinates inside the triangle. exactly zero
        # travel time is checked by test_prt_zero_length_exits.py.
        assert 0.0 <= trk.iloc[2].t < 1e-6

    elif case == "exit":
        # each tracking time is reported in the cell being left,
        # immediately followed by the exit, at the same time
        for t in tracktimes:
            i = trk.index[(trk.ireason == USERTIME) & (trk.t == t)][0]
            nxt = trk.loc[i + 1]
            assert nxt.ireason == EXIT, trk.loc[i - 1 : i + 2]
            assert nxt.t == pytest.approx(t, rel=1e-12)
            assert nxt.icell == trk.loc[i].icell

    elif case == "tstep":
        # each tracking time is reported in the step it ends, before
        # (or with no) time step event at the same time
        for k, t in enumerate(tracktimes, start=1):
            ut = trk[(trk.ireason == USERTIME) & (trk.t == t)].iloc[0]
            assert (ut.kper, ut.kstp) == (1, k)
            ts = trk[(trk.ireason == TIMESTEP) & (trk.t == t)]
            if len(ts):
                assert ts.index[0] > ut.name
                assert (ts.iloc[0].kper, ts.iloc[0].kstp) == (1, k)


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    case, grid = parse(name)
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, grid, case),
        check=lambda t: check_output(t, grid, case),
        targets=targets,
        compare=None,
    )
    test.run()


def test_stay_tstep(function_tmpdir, targets, monkeypatch):
    """Dry cell, STAY method: a tracking time at a period end is reported once."""
    import test_prt_dry as dry

    # end of period 1 (one step)
    monkeypatch.setattr(dry, "user_time", float(dry.perlen[0]))
    name = dry.simname + "_stay"
    idx = dry.cases.index(name)
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: dry.build_models(
            idx, t, newton=True, dry_tracking_method="stay"
        ),
        targets=targets,
        compare=None,
    )
    test.run()
    prt_ws = test.workspace / "prt"
    prtname = get_model_name(name, "prt")
    trk = pd.DataFrame(
        np.fromfile(
            prt_ws / f"{prtname}.trk",
            dtype=get_track_dtype(prt_ws / f"{prtname}.trk.hdr"),
        )
    )
    nparts = trk.irpt.nunique()
    ut = trk[trk.ireason == USERTIME]
    assert len(ut) == nparts
    assert (ut.t == dry.user_time).all()
    assert ut.groupby("irpt").size().eq(1).all()
    assert ((ut.kper == 1) & (ut.kstp == 1)).all()


def test_stay_arrival(function_tmpdir, targets):
    """
    STAY method: a tracking time reported before the particle enters a dry
    cell, partway through a time step, is not reported again.
    """
    mf6 = targets["mf6"]
    name = "tiesarr"
    gwfname, prtname = get_model_name(name, "gwf"), get_model_name(name, "prt")
    gwf_ws, prt_ws = function_tmpdir / "gwf", function_tmpdir / "prt"
    ncol, delr, perlen = 10, 10.0, 1000.0
    dis = dict(nlay=2, nrow=1, ncol=ncol, delr=delr, delc=1.0, top=20.0)
    dis["botm"] = [10.0, 0.0]

    # one row, two layers. water enters the top layer at the left and
    # leaves the bottom layer at the right, so with newton the top layer
    # dries toward the right but still takes inflow from its wet neighbor.
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=gwf_ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[(perlen, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(gwf, **dis)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=1, k33=0.1, save_specific_discharge=True, save_saturation=True
    )
    flopy.mf6.ModflowGwfic(gwf, strt=15.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 15.0], [(1, 0, ncol - 1), 2.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(
        sim,
        linear_acceleration="BICGSTAB",
        outer_maximum=200,
        under_relaxation="DBD",
    )
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]

    tracktimes = [1e-3]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, sim_ws=prt_ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[(perlen, 1, 1.0)])
    prt = flopy.mf6.ModflowPrt(sim, modelname=prtname)
    flopy.mf6.ModflowPrtdis(prt, **dis)
    flopy.mf6.ModflowPrtmip(prt, porosity=0.1)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prtname}_1.prp",
        nreleasepts=1,
        packagedata=[[0, (0, 0, 4), 4.5 * delr, 0.5, 10.37]],
        perioddata={0: [("FIRST",)]},
        dry_tracking_method="stay",
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[f"{prtname}.trk"],
        trackcsv_filerecord=[f"{prtname}.trk.csv"],
        ntracktimes=len(tracktimes),
        tracktimes=[(t,) for t in tracktimes],
    )
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_ws / f"{gwfname}.hds"),
            ("GWFBUDGET", gwf_ws / f"{gwfname}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prtname}.ems")
    sim.register_solution_package(ems, [prt.name])
    sim.write_simulation(silent=True)
    assert sim.run_simulation(silent=True)[0]

    # release just under the water table in the last wet cell of the top
    # layer. the particle leaves it partway through the step, after the
    # tracking time, into the dry cell next door, and stays there.
    trk = read_track(prt_ws, name)
    ex = trk[trk.ireason == EXIT]
    assert len(ex) == 1
    assert tracktimes[0] < ex.iloc[0].t < perlen
    assert (trk.loc[ex.index[0] + 1 :].icell == ex.iloc[0].icell + 1).all()
    check_usertimes_once(trk, tracktimes)

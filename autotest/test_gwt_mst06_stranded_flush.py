"""
Tests the MST stranded mass option for solute that evapotranspiration leaves
in the drained part of a cell, and that recharge flushes from it.

A single cell is held at a saturation of 0.5 by evapotranspiration and
recharge, or by a well, that balance. Evapotranspiration given a mixed
(AUXMIXED) concentration removes water and leaves solute behind; the share in
the drained part of the cell is deposited there, and recharge passing through
the drained pore volume flushes it back as a well-mixed volume.

Cases:
  - steady : evapotranspiration and clean recharge balance for 2,000 d; the
             deposited mass and the concentration reach the steady state in
             which deposition equals flushing, 5,000 g and 66.7 g/m3, for any
             time step. Recharge alone, balanced by a well, then flushes the
             deposit as an exact exponential.
  - et_conc: evapotranspiration at 20 g/m3, balanced by a well injecting clean
             water; only the solute that the mixing of evapotranspiration does
             not carry out of the cell is left behind and deposited.
  - restart: a simulation started from the concentration and the deposit
             (STRANDED_ET FILEOUT, STRANDED_ET_SOLUTE) that another wrote at the
             end of the deposition ends where a single simulation does.
  - errors : STRANDED_ET_SOLUTE without STRANDED_MASS, and a negative deposit,
             are errors (xfail).
  - sorption: a cell with linear sorption drains by half, deposits solute
             under evapotranspiration and recharge, and is refilled to 0.9;
             8,000 g of sorbate is stranded, four fifths of the sorbate and of
             the deposit return, and the cell keeps the solute that the well
             did not remove.
  - restart_all: the same with half the specific yield and first-order
             decay, so that all three reservoirs hold mass, started again
             from the files the first part wrote, ends where a single
             simulation does.
  - ats    : evapotranspiration and recharge with a transport solver that
             rejects the first, long time step, so that ATS repeats it; the
             solute of the cell, in the water and deposited, is conserved.
  - rch_conc: recharge at 10 g/m3, twice the evapotranspiration, balanced by
             a well; the concentration settles at 20 g/m3 and the deposit at
             750 g, so the recharge solute and the flushing are not counted
             twice.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwt_mst06_stranded_cell import (
    assert_budget_closes,
    held_mass,
    net_in,
    strand_budget,
)

delr = delc = 10.0
area = delr * delc
top, botm = 10.0, 0.0
vcell = area * (top - botm)
porosity = 0.3
cinit = 100.0
strt = 5.0
sat = (strt - botm) / (top - botm)
# evapotranspiration at the starting head, which an extinction depth equal to
# the thickness of the cell makes half of the maximum rate
et_max = 0.1
q_et = et_max * area * (strt - (top - (top - botm))) / (top - botm)
vdrained = porosity * (1.0 - sat) * vcell


def get_model(
    ws,
    periods,
    c0=cinit,
    strand=True,
    mst_kwargs=None,
    h0=strt,
    sy=porosity,
    rcha=False,
    ats=None,
    gwt_ims=None,
    newton=False,
):
    """Single cell with evapotranspiration, recharge, and a well.

    Each period is (length, number of time steps, maximum evapotranspiration
    rate, concentration of evapotranspiration, recharge rate, well rate), and
    optionally the concentration of recharge, zero if not given. With rcha,
    recharge is given as an array (RCHA); ats is the ATS period data,
    gwt_ims overrides the settings of the transport solver, and newton selects
    the Newton formulation, which keeps a cell drained from its top by a well
    active.
    """
    sim = flopy.mf6.MFSimulation(sim_name="flush", sim_ws=ws, exe_name="mf6")
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(p[0], p[1], 1.0) for p in periods],
    )
    if ats is not None:
        tdis.ats.initialize(maxats=len(ats), perioddata=ats, filename="flush.ats")

    def add_ims(fname, **settings):
        kwargs = {
            "complexity": "moderate",
            "outer_dvclose": 1e-10,
            "inner_dvclose": 1e-12,
            "linear_acceleration": "bicgstab",
        }
        kwargs.update(settings)
        return flopy.mf6.ModflowIms(sim, filename=fname, **kwargs)

    def c_rch(p):
        return p[6] if len(p) > 6 else 0.0

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="gwf",
        save_flows=True,
        newtonoptions="NEWTON" if newton else None,
    )
    sim.register_ims_package(add_ims("gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=h0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0, save_saturation=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=sy, transient={0: True})
    cell = (0, 0, 0)
    flopy.mf6.ModflowGwfevt(
        gwf,
        stress_period_data={
            k: [[cell, top, p[2], top - botm, p[3]]] for k, p in enumerate(periods)
        },
        auxiliary="CONCENTRATION",
        pname="EVT-1",
    )
    if rcha:
        flopy.mf6.ModflowGwfrcha(
            gwf,
            recharge={k: p[4] for k, p in enumerate(periods)},
            auxiliary="CONCENTRATION",
            aux={k: [c_rch(p)] for k, p in enumerate(periods)},
            pname="RCH-1",
        )
    else:
        flopy.mf6.ModflowGwfrch(
            gwf,
            stress_period_data={
                k: [[cell, p[4], c_rch(p)]] for k, p in enumerate(periods)
            },
            auxiliary="CONCENTRATION",
            pname="RCH-1",
        )
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={k: [[cell, p[5], 0.0]] for k, p in enumerate(periods)},
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    sim.register_ims_package(add_ims("gwt.ims", **(gwt_ims or {})), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=c0)
    kwargs = {"porosity": porosity, "save_flows": True}
    if strand:
        kwargs.update(stranded_mass=True, stranded_filerecord=[("gwt.strand.bin",)])
    if mst_kwargs is not None:
        kwargs.update(mst_kwargs)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[
            ["EVT-1", "AUXMIXED", "CONCENTRATION"],
            ["RCH-1", "AUX", "CONCENTRATION"],
            ["WEL-1", "AUX", "CONCENTRATION"],
        ],
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def series(ws, fname, text):
    f = flopy.utils.HeadFile(ws / fname, text=text, precision="double")
    times = np.array(f.get_times())
    return times, np.array([f.get_data(totim=t).ravel()[0] for t in times])


def run_framework(function_tmpdir, targets, build, check):
    TestFramework(
        name="flush",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


@pytest.mark.parametrize("rcha", [False, True], ids=["rch", "rcha"])
@pytest.mark.parametrize("nstp", [20, 200])
def test_steady(function_tmpdir, targets, nstp, rcha):
    # clean recharge balances evapotranspiration, then recharge alone is
    # balanced by a well that removes water at the concentration of the cell
    q_rch = q_et
    periods = [
        (2000.0, nstp, et_max, 0.0, q_rch / area, 0.0),
        (300.0, 30, 0.0, 0.0, q_rch / area, -q_rch),
    ]

    def build(test):
        return get_model(test.workspace, periods, rcha=rcha)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, head = series(ws, "gwf.hds", "HEAD")
        assert np.allclose(head, strt), f"head {head}"
        times, conc = series(ws, "gwt.ucn", "CONCENTRATION")
        _, stranded = series(ws, "gwt.strand.bin", "STRANDED")
        # the mass of the cell, in the water and deposited, does not change
        # while clean recharge balances evapotranspiration; at steady state
        # deposition equals flushing, so the deposit is (1 - S) Q_et C V_d /
        # Q_rch, which is 75 C here
        held_ss = (1.0 - sat) * q_et * vdrained / q_rch
        mobile = porosity * sat * vcell
        c_ss = cinit * mobile / (mobile + held_ss)
        end1 = np.isclose(times, periods[0][0])
        assert np.isclose(conc[end1][0], c_ss, rtol=1e-6), (
            f"concentration {conc[end1][0]} expected {c_ss}"
        )
        assert np.isclose(stranded[end1][0], held_ss * c_ss, rtol=1e-6), (
            f"deposit {stranded[end1][0]} expected {held_ss * c_ss}"
        )
        # recharge alone flushes the deposit as an exact exponential
        flush = times > periods[0][0]
        k = q_rch / vdrained
        expected = stranded[end1][0] * np.exp(-k * (times[flush] - periods[0][0]))
        assert np.allclose(stranded[flush], expected, rtol=1e-8), (
            f"deposit {stranded[flush]} expected {expected}"
        )
        # the deposit is reported apart from drainage
        strand, _ = strand_budget(ws)
        assert np.allclose(held_mass(strand, "AQUEOUS"), 0.0, atol=1e-6)
        assert net_in(strand, "ET-DEPOSIT").to_numpy()[-1] > 0.0
        assert net_in(strand, "RECHARGE-FLUSH").to_numpy()[-1] < 0.0

    run_framework(function_tmpdir, targets, build, check)


def test_et_conc(function_tmpdir, targets):
    c_et = 20.0
    periods = [(100.0, 10, et_max, c_et, 0.0, q_et)]

    def build(test):
        return get_model(test.workspace, periods)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, conc = series(ws, "gwt.ucn", "CONCENTRATION")
        assert conc[-1] > c_et, f"concentration {conc}"
        # what the mixing of evapotranspiration does not carry out of the cell
        # is left behind, and the share in the drained part is deposited
        lb = flopy.utils.Mf6ListBudget(ws / "gwt.lst", budgetkey="MASS BUDGET FOR MST")
        deposit = net_in(
            lb.get_dataframes(diff=False, start_datetime=None)[0], "ET-DEPOSIT"
        ).to_numpy()
        lbm = flopy.utils.Mf6ListBudget(
            ws / "gwt.lst", budgetkey="MASS BUDGET FOR ENTIRE MODEL"
        )
        incm = lbm.get_dataframes(diff=False, start_datetime=None)[0]
        evt_out = incm["EVT_OUT"].to_numpy()
        expected = (1.0 - sat) * (q_et * conc - evt_out)
        assert np.allclose(deposit, expected, rtol=1e-5), (
            f"deposit {deposit} expected {expected}"
        )
        assert deposit.min() > 0.0, f"deposit {deposit}"

    run_framework(function_tmpdir, targets, build, check)


def test_restart(function_tmpdir, targets):
    exe = str(targets["mf6"])
    deposit = (500.0, 50, et_max, 0.0, q_et / area, 0.0)
    flush = (300.0, 30, 0.0, 0.0, q_et / area, -q_et)
    dep_out = {"stranded_et_filerecord": [("gwt.et.bin",)]}

    def run(ws, periods, **kwargs):
        sim = get_model(ws, periods, **kwargs)
        sim.exe_name = exe
        sim.write_simulation(silent=True)
        success, buff = sim.run_simulation(silent=True)
        assert success, f"simulation in {ws} failed\n{buff}"

    # one simulation over both halves, and the first half alone
    run(function_tmpdir / "full", [deposit, flush], mst_kwargs=dep_out)
    run(function_tmpdir / "first", [deposit], mst_kwargs=dep_out)
    _, c1 = series(function_tmpdir / "first", "gwt.ucn", "CONCENTRATION")
    _, d1 = series(function_tmpdir / "first", "gwt.et.bin", "STRAND-ET")
    assert d1[-1] > 0.0, f"no solute was deposited {d1}"

    # the second half, started from what the first half wrote
    run(
        function_tmpdir / "second",
        [flush],
        c0=c1[-1],
        mst_kwargs={"stranded_et_solute": d1[-1]},
    )
    for fname, text in (
        ("gwt.ucn", "CONCENTRATION"),
        ("gwt.strand.bin", "STRANDED"),
    ):
        _, full = series(function_tmpdir / "full", fname, text)
        _, second = series(function_tmpdir / "second", fname, text)
        assert np.allclose(second[-1], full[-1], rtol=1e-8), (
            f"{text} {second[-1]} after a restart, {full[-1]} without"
        )
    # the deposit file holds the deposit alone, which is all that is stranded
    _, dfull = series(function_tmpdir / "full", "gwt.et.bin", "STRAND-ET")
    _, sfull = series(function_tmpdir / "full", "gwt.strand.bin", "STRANDED")
    assert np.allclose(dfull, sfull), "the deposit is not the stranded mass"


deposit_errors = {
    "nostrand": (
        False,
        10.0,
        "STRANDED_AQUEOUS, STRANDED_SORBED, or STRANDED_ET_SOLUTE was specified "
        "but the STRANDED_MASS keyword was not specified",
    ),
    "negative": (
        True,
        -10.0,
        "STRANDED_AQUEOUS, STRANDED_SORBED, and STRANDED_ET_SOLUTE cannot be negative",
    ),
}


@pytest.mark.parametrize("case", list(deposit_errors))
def test_errors(function_tmpdir, targets, case):
    strand, value, message = deposit_errors[case]
    periods = [(10.0, 1, et_max, 0.0, 0.0, q_et)]

    def build(test):
        return get_model(
            test.workspace,
            periods,
            strand=strand,
            mst_kwargs={"stranded_et_solute": value},
        )

    def check(test):
        text = " ".join((test.workspace / "mfsim.lst").read_text().split())
        assert message in text, f"expected {message!r}"

    TestFramework(
        name="flush",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        xfail=True,
    ).run()


rhob, kd = 1600.0, 1.0e-4
retard = porosity + rhob * kd


def test_sorption(function_tmpdir, targets):
    # the well drains the saturated cell to half in a day, evapotranspiration
    # and recharge then deposit and flush solute at that saturation, and
    # clean water injected for a day refills the cell to 0.9
    q_drain = porosity * 0.5 * vcell
    q_refill = porosity * 0.4 * vcell
    periods = [
        (1.0, 1, 0.0, 0.0, 0.0, -q_drain),
        (500.0, 50, et_max, 0.0, q_et / area, 0.0),
        (1.0, 1, 0.0, 0.0, 0.0, q_refill),
    ]

    def build(test):
        return get_model(
            test.workspace,
            periods,
            h0=top,
            newton=True,
            mst_kwargs={"sorption": "linear", "bulk_density": rhob, "distcoef": kd},
        )

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, conc = series(ws, "gwt.ucn", "CONCENTRATION")
        _, stranded = series(ws, "gwt.strand.bin", "STRANDED")
        strand, _ = strand_budget(ws)
        srb = held_mass(strand, "SORBED")
        et = held_mass(strand, "ET")
        # the specific yield equals the porosity, so only the sorbate of the
        # drained half of the cell, at 100 g/m3, is stranded by the well
        expected = 0.5 * vcell * rhob * kd * cinit
        assert np.isclose(srb[0], expected, rtol=1e-6), f"sorbed {srb[0]}"
        # the budget is printed every time step, so the second last entry is
        # the end of the deposition, before the refill
        assert np.isclose(srb[-2], expected, rtol=1e-6), f"sorbed {srb[-2]}"
        assert et[-2] > 0.0, f"no solute was deposited {et}"
        # the refill resaturates 0.4 of the 0.5 that drained, so four fifths
        # of the sorbate and of the deposit return
        expected = 0.2 * (srb[-2] + et[-2])
        assert np.isclose(stranded[-1], expected, rtol=1e-5), (
            f"stranded {stranded[-1]} expected {expected}"
        )
        # evapotranspiration and recharge carry no solute, so the cell, in
        # the water, on the solids, and held, keeps all but the solute the
        # well removed
        total = retard * 0.9 * vcell * conc[-1] + stranded[-1]
        expected = retard * vcell * cinit - q_drain * cinit
        assert np.isclose(total, expected, rtol=1e-6), (
            f"solute of the cell {total} expected {expected}"
        )

    run_framework(function_tmpdir, targets, build, check)


def test_restart_all(function_tmpdir, targets):
    exe = str(targets["mf6"])
    sy = 0.5 * porosity
    q_drain = sy * 0.5 * vcell
    drain = (1.0, 1, 0.0, 0.0, 0.0, -q_drain)
    deposit = (300.0, 30, et_max, 0.0, q_et / area, 0.0)
    # refill to a saturation of 0.75, so that half of each reservoir returns
    refill = (1.0, 1, 0.0, 0.0, 0.0, 0.5 * q_drain)
    mst = {
        "sorption": "linear",
        "bulk_density": rhob,
        "distcoef": kd,
        "first_order_decay": True,
        "decay": 0.002,
        "decay_sorbed": 0.001,
    }
    out = {
        "stranded_aq_filerecord": [("gwt.aq.bin",)],
        "stranded_srb_filerecord": [("gwt.srb.bin",)],
        "stranded_et_filerecord": [("gwt.et.bin",)],
    }

    def run(ws, periods, **kwargs):
        sim = get_model(ws, periods, sy=sy, newton=True, **kwargs)
        sim.exe_name = exe
        sim.write_simulation(silent=True)
        success, buff = sim.run_simulation(silent=True)
        assert success, f"simulation in {ws} failed\n{buff}"

    run(function_tmpdir / "full", [drain, deposit, refill], h0=top, mst_kwargs=mst)
    run(function_tmpdir / "first", [drain, deposit], h0=top, mst_kwargs=mst | out)
    first = function_tmpdir / "first"
    _, c1 = series(first, "gwt.ucn", "CONCENTRATION")
    held = {
        name: series(first, fname, text)[1][-1]
        for name, fname, text in (
            ("stranded_aqueous", "gwt.aq.bin", "STRAND-AQUEOUS"),
            ("stranded_sorbed", "gwt.srb.bin", "STRAND-SORBED"),
            ("stranded_et_solute", "gwt.et.bin", "STRAND-ET"),
        )
    }
    assert all(v > 0.0 for v in held.values()), f"held at the restart {held}"

    run(
        function_tmpdir / "second",
        [refill],
        c0=c1[-1],
        h0=strt,
        mst_kwargs=mst | held,
    )
    assert_budget_closes(function_tmpdir / "full")
    for fname, text in (
        ("gwt.ucn", "CONCENTRATION"),
        ("gwt.strand.bin", "STRANDED"),
    ):
        _, full = series(function_tmpdir / "full", fname, text)
        _, second = series(function_tmpdir / "second", fname, text)
        assert np.isclose(second[-1], full[-1], rtol=1e-6), (
            f"{text} {second[-1]} after a restart, {full[-1]} without"
        )


def test_ats(function_tmpdir, targets):
    # a transport solution allowed one outer iteration converges only if the
    # concentration changes by less than 20 g/m3, so ATS repeats the first,
    # 50 d time step with a fifth of its length; the steps are limited to
    # between 1 and 10 d after that, and a short last period keeps a repeated
    # step from being the last of the simulation
    q_rch = q_et
    periods = [
        (200.0, 1, et_max, 0.0, q_rch / area, 0.0),
        (1.0, 1, et_max, 0.0, q_rch / area, 0.0),
    ]

    def build(test):
        return get_model(
            test.workspace,
            periods,
            ats=[(0, 50.0, 1.0, 50.0, 1.0, 5.0)],
            gwt_ims={"outer_maximum": 1, "outer_dvclose": 20.0},
        )

    def check(test):
        ws = test.workspace
        listing = (ws / "mfsim.lst").read_text()
        assert "Failed solution for step" in listing, "no time step was repeated"
        times, conc = series(ws, "gwt.ucn", "CONCENTRATION")
        _, stranded = series(ws, "gwt.strand.bin", "STRANDED")
        assert len(times) < 50, f"{len(times)} time steps"
        # the solute of the cell, in the water and deposited, is conserved
        # through the repeated time step, and solute was deposited
        mobile = porosity * sat * vcell
        total = mobile * conc + stranded
        assert np.allclose(total, mobile * cinit, rtol=1e-9), f"total {total}"
        assert stranded[-1] > 0.0, f"no solute was deposited {stranded}"

    run_framework(function_tmpdir, targets, build, check)


def test_rch_conc(function_tmpdir, targets):
    # recharge at twice the evapotranspiration and 10 g/m3, and a well that
    # removes the rest of the water at the concentration of the cell
    c_rch = 10.0
    q_rch = 2.0 * q_et
    q_wel = -(q_rch - q_et)
    periods = [(2000.0, 200, et_max, 0.0, q_rch / area, q_wel, c_rch)]

    def build(test):
        return get_model(test.workspace, periods)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, conc = series(ws, "gwt.ucn", "CONCENTRATION")
        _, stranded = series(ws, "gwt.strand.bin", "STRANDED")
        # at steady state the deposit and the flushing balance, so the water
        # balance alone sets the concentration, Q_rch C_rch = Q_wel C
        c_ss = q_rch * c_rch / -q_wel
        held_ss = (1.0 - sat) * q_et * c_ss * vdrained / q_rch
        assert np.isclose(conc[-1], c_ss, rtol=1e-6), f"{conc[-1]} expected {c_ss}"
        assert np.isclose(stranded[-1], held_ss, rtol=1e-6), (
            f"deposit {stranded[-1]} expected {held_ss}"
        )

    run_framework(function_tmpdir, targets, build, check)

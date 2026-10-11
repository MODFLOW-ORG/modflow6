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


def get_model(ws, periods):
    """Single cell with evapotranspiration, recharge, and a well.

    Each period is (length, number of time steps, maximum evapotranspiration
    rate, concentration of evapotranspiration, recharge rate, well rate).
    """
    sim = flopy.mf6.MFSimulation(sim_name="flush", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(p[0], p[1], 1.0) for p in periods],
    )

    def add_ims(fname):
        return flopy.mf6.ModflowIms(
            sim,
            complexity="moderate",
            outer_dvclose=1e-10,
            inner_dvclose=1e-12,
            linear_acceleration="bicgstab",
            filename=fname,
        )

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    sim.register_ims_package(add_ims("gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0, save_saturation=True)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=porosity, transient={0: True})
    cell = (0, 0, 0)
    flopy.mf6.ModflowGwfevt(
        gwf,
        stress_period_data={
            k: [[cell, top, p[2], top - botm, p[3]]] for k, p in enumerate(periods)
        },
        auxiliary="CONCENTRATION",
        pname="EVT-1",
    )
    flopy.mf6.ModflowGwfrch(
        gwf,
        stress_period_data={k: [[cell, p[4], 0.0]] for k, p in enumerate(periods)},
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
    sim.register_ims_package(add_ims("gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        save_flows=True,
        stranded_mass=True,
        stranded_filerecord=[("gwt.strand.bin",)],
    )
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


@pytest.mark.parametrize("nstp", [20, 200])
def test_steady(function_tmpdir, targets, nstp):
    # clean recharge balances evapotranspiration, then recharge alone is
    # balanced by a well that removes water at the concentration of the cell
    q_rch = q_et
    periods = [
        (2000.0, nstp, et_max, 0.0, q_rch / area, 0.0),
        (300.0, 30, 0.0, 0.0, q_rch / area, -q_rch),
    ]

    def build(test):
        return get_model(test.workspace, periods)

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

"""
Tests the MST stranded mass option against analytical solutions for a single
cell that a well drains and rewets.

The cell is 10 m on each side and starts at 100 g/m3, the example used by
Bedekar and others (written communication) to compare treatments of sorbed mass in a
cell whose water table moves. With the option, mass held out of the mobile
domain is returned when the cell rewets, so the solutions do not depend on the
length of the time step except where the concentration changes while the cell
drains.

Cases:
  - paper          : pump 25 m3/d for 2 d and inject clean water for 2 d, with
                     specific yield equal to the porosity; 50 g/m3 at the end
                     without sorption, and 7,500 g stranded, 20,000 g in the
                     cell, and 80 g/m3 at the end with sorption.
  - retained       : the same with half the specific yield and no sorption; the
                     water retained against drainage strands 2,500 g, and the
                     cell ends at 75 g/m3.
  - dilution       : the cell drains while clean water is injected; the
                     concentration and stranded mass converge to the analytical
                     solution as the time step is refined.
  - confined_start : the head starts above the top of the cell, with specific
                     storage; the mass stranded follows the STO-SY flow, and
                     the reservoir empties when the head returns above the top.
  - rise_above     : the cell starts partly saturated, drains, and rises above
                     its starting level into material that never drained; no
                     mass is created.
  - decay          : sorbed mass stranded by pumping decays at the sorbed
                     rate while the cell stays drained, as an exact
                     exponential.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

top, botm = 10.0, 0.0
delr = delc = 10.0
area = delr * delc
vcell = area * (top - botm)
porosity = 0.1
rhob = 1500.0
kd = 1.0e-4
rhobkd = rhob * kd
cinit = 100.0
# pump, then inject clean water, at 25 m3/d; each entry is (length, wells)
paper = [(2.0, [(-25.0, 0.0)]), (2.0, [(25.0, 0.0)])]


def get_model(
    ws,
    periods,
    nstp,
    sy,
    sorb,
    strt=top,
    ss=0.0,
    strand=True,
    mst_kwargs=None,
    transport=True,
):
    """Single cell drained and rewetted by wells, with a GWF-GWT exchange.

    Without transport, only the flow model is built.
    """
    sim = flopy.mf6.MFSimulation(sim_name="cell", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(perlen, nstp, 1.0) for perlen, _ in periods],
    )

    # the moderate settings damp the change in storage coefficient at the top
    # of the cell, where the cell starts
    def add_ims(fname):
        return flopy.mf6.ModflowIms(
            sim,
            complexity="moderate",
            outer_dvclose=1e-10,
            inner_dvclose=1e-12,
            linear_acceleration="bicgstab",
            filename=fname,
        )

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", save_flows=True, newtonoptions="NEWTON"
    )
    sim.register_ims_package(add_ims("gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    # specific discharge and saturation are read by FMI when transport runs
    # in a separate simulation
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=10.0,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=ss, sy=sy, transient={0: True})
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={
            k: [[(0, 0, 0), q, c] for q, c in wells]
            for k, (_, wells) in enumerate(periods)
        },
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    if not transport:
        return sim

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    sim.register_ims_package(add_ims("gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    kwargs = {"porosity": porosity, "save_flows": True}
    if sorb:
        kwargs.update(sorption="linear", bulk_density=rhob, distcoef=kd)
    if strand:
        kwargs.update(stranded_mass=True, stranded_filerecord=[("gwt.strand.bin",)])
    if mst_kwargs is not None:
        kwargs.update(mst_kwargs)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def results(ws, sorb):
    """Times, heads, concentrations, stranded mass, and mass in the cell."""
    hf = flopy.utils.HeadFile(os.path.join(ws, "gwf.hds"))
    times = np.array(hf.get_times())
    head = np.array([hf.get_data(totim=t).flatten()[0] for t in times])
    sat = np.clip((head - botm) / (top - botm), 0.0, 1.0)
    cf = flopy.utils.HeadFile(os.path.join(ws, "gwt.ucn"), text="CONCENTRATION")
    conc = np.array([cf.get_data(totim=t).flatten()[0] for t in times])
    sf = flopy.utils.HeadFile(os.path.join(ws, "gwt.strand.bin"), text="STRANDED")
    stranded = np.array([sf.get_data(totim=t).flatten()[0] for t in times])
    retard = porosity + (rhobkd if sorb else 0.0)
    mobile = retard * sat * vcell * conc
    return times, head, sat, conc, stranded, mobile + stranded


def assert_budget_closes(ws):
    """Cumulative percent discrepancies of the transport model budget are zero."""
    disc = []
    model_budget = False
    with open(os.path.join(ws, "gwt.lst")) as f:
        for line in f:
            if "BUDGET FOR ENTIRE MODEL" in line:
                model_budget = True
            elif "PERCENT DISCREPANCY" in line and model_budget:
                disc.append(
                    float(line.replace("PERCENT DISCREPANCY =", " ").split()[0])
                )
                model_budget = False
    assert disc, "no model budget found in gwt.lst"
    assert np.allclose(disc, 0.0, atol=1e-2), f"model budget does not close {disc}"


def run_framework(function_tmpdir, targets, build, check):
    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


@pytest.mark.parametrize("nstp", [1, 10])
@pytest.mark.parametrize("sorb", [False, True])
def test_paper(function_tmpdir, targets, sorb, nstp):
    def build(test):
        return get_model(test.workspace, paper, nstp, porosity, sorb)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, _, sat, conc, stranded, mass = results(ws, sorb)
        mid = nstp - 1
        # pumping removes 5,000 g of dissolved mass at 100 g/m3; with sorption
        # the 7,500 g sorbed on the half of the cell that drains is stranded
        assert np.isclose(sat[mid], 0.5), f"saturation after pumping {sat[mid]}"
        assert np.isclose(conc[mid], cinit), f"concentration {conc[mid]}"
        assert np.isclose(stranded[mid], 7500.0 if sorb else 0.0, atol=1e-6)
        assert np.isclose(mass[mid], 20000.0 if sorb else 5000.0)
        # the clean water dilutes 5,000 g in 100 m3, or 20,000 g in the water
        # and sorbate of the full cell, and returns all of the stranded mass
        assert np.isclose(sat[-1], 1.0), f"saturation after injection {sat[-1]}"
        assert np.isclose(conc[-1], 80.0 if sorb else 50.0), f"{conc[-1]}"
        assert np.isclose(stranded[-1], 0.0, atol=1e-6), f"{stranded[-1]}"
        assert np.isclose(mass[-1], 20000.0 if sorb else 5000.0), f"{mass[-1]}"

    run_framework(function_tmpdir, targets, build, check)


@pytest.mark.parametrize("nstp", [1, 10])
def test_retained(function_tmpdir, targets, nstp):
    # half of the drained pore space is retained against drainage
    periods = [(1.0, [(-25.0, 0.0)]), (1.0, [(25.0, 0.0)])]

    def build(test):
        return get_model(test.workspace, periods, nstp, 0.5 * porosity, False)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, _, sat, conc, stranded, mass = results(ws, False)
        mid = nstp - 1
        # 25 m3 drains from 50 m3 of pore space, so 25 m3 at 100 g/m3 stays
        assert np.isclose(sat[mid], 0.5), f"saturation after pumping {sat[mid]}"
        assert np.isclose(conc[mid], cinit), f"concentration {conc[mid]}"
        assert np.isclose(stranded[mid], 2500.0), f"stranded {stranded[mid]}"
        # 7,500 g remain after 2,500 g are pumped out, in 100 m3 of water
        assert np.isclose(conc[-1], 75.0), f"concentration {conc[-1]}"
        assert np.isclose(stranded[-1], 0.0, atol=1e-6), f"{stranded[-1]}"
        assert np.isclose(mass[-1], 7500.0), f"mass {mass[-1]}"

    run_framework(function_tmpdir, targets, build, check)


def test_dilution(function_tmpdir, targets):
    """Stranded mass converges to the analytical solution with the time step.

    A well removes 50 m3/d while another injects 25 m3/d of clean water, so the
    cell drains while it is diluted. For a well-mixed cell the concentration is
    C0 (S/S0)**p, with p = Qin theta / (R Q), and the stranded mass is
    rhob Kd V C0 S0 (1 - (S/S0)**(p + 1)) / (p + 1).
    """
    exe = str(targets["mf6"])
    s0 = 0.9
    qin, qnet = 25.0, 25.0
    periods = [(2.0, [(-50.0, 0.0), (qin, 0.0)])]
    p = qin * porosity / ((porosity + rhobkd) * qnet)
    ratio = (s0 - qnet * 2.0 / (porosity * vcell)) / s0
    cexact = cinit * ratio**p
    sexact = rhobkd * vcell * cinit * s0 * (1.0 - ratio ** (p + 1.0)) / (p + 1.0)

    errors = []
    for nstp in (4, 16, 64):
        ws = function_tmpdir / f"nstp{nstp}"
        sim = get_model(ws, periods, nstp, porosity, True, strt=s0 * top)
        sim.exe_name = exe
        sim.write_simulation(silent=True)
        success, buff = sim.run_simulation(silent=True)
        assert success, f"nstp={nstp} failed\n{buff}"
        assert_budget_closes(ws)
        _, _, _, conc, stranded, _ = results(ws, True)
        errors.append((abs(conc[-1] / cexact - 1.0), abs(stranded[-1] / sexact - 1.0)))

    # backward Euler is first order, so each fourfold refinement cuts the
    # error by about four
    errors = np.array(errors)
    assert np.all(errors[1:] < 0.3 * errors[:-1]), f"not converging {errors}"
    assert np.all(errors[-1] < 5e-3), f"error at 64 steps {errors[-1]}"


@pytest.mark.parametrize("nstp", [1, 4])
def test_confined_start(function_tmpdir, targets, nstp):
    # the head starts 2 m above the top of the cell, falls below it, and
    # returns; specific storage releases water above and below the top
    sy = 0.5 * porosity
    periods = [(1.0, [(-25.0, 0.0)]), (1.0, [(27.0, 0.0)])]

    def build(test):
        return get_model(
            test.workspace, periods, nstp, sy, True, strt=top + 2.0, ss=1.0e-3
        )

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, head, sat, conc, stranded, _ = results(ws, True)
        # only the water released by specific yield drains pore space; the
        # rest of the change in the mobile water volume is retained
        cbc = flopy.utils.CellBudgetFile(
            os.path.join(ws, "gwf.cbc"), precision="double"
        )
        stosy = np.array([np.ravel(r)[0] for r in cbc.get_data(text="STO-SY")])
        delt = 1.0 / nstp
        sat_old = np.concatenate(([1.0], sat[:-1]))
        drain = np.clip(sat_old - sat, 0.0, None)
        retained = np.clip(drain * porosity * vcell - stosy * delt, 0.0, None)
        expected = np.cumsum(retained * conc + drain * rhobkd * vcell * conc)
        assert sat[nstp - 1] < 1.0, "the cell did not drain"
        # the saturation computed from the heads differs slightly from the one
        # the flow model uses with the Newton formulation
        assert np.allclose(stranded[:nstp], expected[:nstp], rtol=1e-4), (
            f"stranded {stranded[:nstp]} expected {expected[:nstp]}"
        )
        assert head[-1] > top, f"the head did not return above the top {head[-1]}"
        assert np.isclose(stranded[-1], 0.0, atol=1e-6), f"{stranded[-1]}"

    run_framework(function_tmpdir, targets, build, check)


def test_rise_above(function_tmpdir, targets):
    # the cell starts half saturated, drains to 0.3, and rises to 0.9
    nstp = 4
    periods = [(1.0, [(-20.0, 0.0)]), (1.0, [(60.0, 0.0)])]

    def build(test):
        return get_model(test.workspace, periods, nstp, porosity, True, strt=5.0)

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        _, _, sat, conc, stranded, mass = results(ws, True)
        assert np.isclose(sat[-1], 0.9), f"saturation at the end {sat[-1]}"
        assert stranded.max() > 0.0, "no mass was stranded"
        assert np.isclose(stranded[-1], 0.0, atol=1e-6), f"{stranded[-1]}"
        # mass leaves only through the pumping well, at the concentration of
        # the cell; the clean water and the material that never drained add
        # none
        pumped = np.sum(20.0 * conc[:nstp] / nstp)
        initial = (porosity + rhobkd) * 0.5 * vcell * cinit
        assert np.isclose(mass[-1], initial - pumped, rtol=1e-6), (
            f"mass {mass[-1]} expected {initial - pumped}"
        )

    run_framework(function_tmpdir, targets, build, check)


def test_decay(function_tmpdir, targets):
    # pump for 2 d, then leave the cell drained for 10 d
    lam, lam_srb = 0.05, 0.02
    periods = [(2.0, [(-25.0, 0.0)]), (10.0, [(0.0, 0.0)])]

    def build(test):
        return get_model(
            test.workspace,
            periods,
            10,
            porosity,
            True,
            mst_kwargs={
                "first_order_decay": True,
                "decay": lam,
                "decay_sorbed": lam_srb,
            },
        )

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        times, _, _, _, stranded, _ = results(ws, True)
        # the specific yield equals the porosity, so only sorbed mass is
        # stranded, and it decays at the sorbed rate once the cell stops
        # draining
        held = times > 2.0
        expected = stranded[~held][-1] * np.exp(-lam_srb * (times[held] - 2.0))
        assert stranded[~held][-1] > 0.0, "no mass was stranded"
        assert np.allclose(stranded[held], expected, rtol=1e-10), (
            f"stranded {stranded[held]} expected {expected}"
        )

    run_framework(function_tmpdir, targets, build, check)

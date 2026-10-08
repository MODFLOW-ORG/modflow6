"""
Tests the initial stranded mass of the MST Package (STRANDED_AQUEOUS and
STRANDED_SORBED), which holds solute in the part of a cell that is drained at
the start of the simulation and returns it as that part resaturates.

A single cell, 10 m on each side and 10 m thick, starts half saturated with no
solute in the water. The specific yield equals the porosity, so no water is
retained against drainage and only the initial stranded mass is held.

Cases:
  - refill   : clean water injected at 25 m3/d raises the saturation from 0.5
               to 0.95; each time step returns the initial mass in proportion
               to the rise over the 0.5 drained at the start, so 90 percent
               is returned and the cell ends at 9 g/m3.
  - decay    : the cell stays drained with aqueous and sorbed stranded mass,
               which decay at their own rates as exact exponentials.
  - fmi      : the refill case with flows read from a budget file, after a
               first stress period in which the water table does not move,
               gives the same result.
  - restart  : a simulation started from the concentration and the aqueous
               stranded mass (STRANDED_AQ FILEOUT) that another wrote halfway
               through ends where a single simulation over both halves does.
  - errors   : the arrays without STRANDED_MASS, STRANDED_SORBED without
               sorption, a negative mass, and stranded mass in a cell that
               starts saturated are errors (xfail).
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

top, botm = 10.0, 0.0
delr = delc = 10.0
vcell = delr * delc * (top - botm)
porosity = 0.1
rhob, kd = 1500.0, 1.0e-4
strt0 = 5.0
mass0 = 1000.0
qin = 25.0
# saturation rises by qin / (sy * area * thickness) per day
dsdt = qin / (porosity * vcell)


def add_flow(sim, strt, periods):
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
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=10.0,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=porosity, transient={0: True})
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={
            k: [[(0, 0, 0), q, 0.0]] for k, (_, q) in enumerate(periods)
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
    return add_ims


def add_transport(sim, add_ims, cinit, mst_kwargs, fmi=None):
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    sim.register_ims_package(add_ims("gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt, nrow=1, ncol=1, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    kwargs = {
        "porosity": porosity,
        "save_flows": True,
        "stranded_mass": True,
        "stranded_filerecord": [("gwt.strand.bin",)],
        "stranded_aq_filerecord": [("gwt.strandaq.bin",)],
        "stranded_srb_filerecord": [("gwt.strandsrb.bin",)],
    }
    kwargs.update(mst_kwargs)
    flopy.mf6.ModflowGwtmst(gwt, **kwargs)
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    if fmi is not None:
        flopy.mf6.ModflowGwtfmi(gwt, packagedata=fmi)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )


def get_model(ws, periods, nstp, mst_kwargs, strt=strt0, cinit=0.0):
    """Single cell with flow and transport solved together."""
    sim = flopy.mf6.MFSimulation(sim_name="cell", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(perlen, nstp, 1.0) for perlen, _ in periods],
    )
    add_ims = add_flow(sim, strt, periods)
    add_transport(sim, add_ims, cinit, mst_kwargs)
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return sim


def read_array(ws, fname, text):
    f = flopy.utils.HeadFile(os.path.join(ws, fname), text=text, precision="double")
    times = np.array(f.get_times())
    return times, np.array([f.get_data(totim=t).flatten()[0] for t in times])


def assert_budget_closes(ws):
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


def check_refill(ws, t0=0.0):
    """Check the refill that starts at time t0."""
    times, conc = read_array(ws, "gwt.ucn", "CONCENTRATION")
    _, held = read_array(ws, "gwt.strandaq.bin", "STRAND-AQUEOUS")
    _, total = read_array(ws, "gwt.strand.bin", "STRANDED")
    keep = times > t0
    times, conc, held, total = times[keep] - t0, conc[keep], held[keep], total[keep]
    # the drained part at the start, 0.5, returns as the saturation rises
    sat = strt0 / (top - botm) + dsdt * times
    expected = mass0 * (1.0 - (sat - 0.5) / 0.5)
    assert np.allclose(held, expected, rtol=1e-6), f"held {held} expected {expected}"
    assert np.allclose(total, held), "STRANDED is not the sum of the reservoirs"
    # the returned mass is the only solute in the water
    mobile = porosity * sat * vcell * conc
    assert np.allclose(mobile + held, mass0, rtol=1e-6), f"mass {mobile + held}"
    assert np.isclose(conc[-1], 0.9 * mass0 / (porosity * 0.95 * vcell), rtol=1e-6)


# inject clean water for 1.8 d, which raises the saturation to 0.95
refill = [(1.8, qin)]


def test_refill(function_tmpdir, targets):
    def build(test):
        return get_model(test.workspace, refill, 9, {"stranded_aqueous": mass0})

    def check(test):
        assert_budget_closes(test.workspace)
        check_refill(test.workspace)

    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


def test_decay(function_tmpdir, targets):
    lam_aq, lam_srb = 0.1, 0.02
    msrb0 = 400.0

    def build(test):
        return get_model(
            test.workspace,
            [(10.0, 0.0)],
            10,
            {
                "sorption": "linear",
                "bulk_density": rhob,
                "distcoef": kd,
                "first_order_decay": True,
                "decay": lam_aq,
                "decay_sorbed": lam_srb,
                "stranded_aqueous": mass0,
                "stranded_sorbed": msrb0,
            },
        )

    def check(test):
        ws = test.workspace
        assert_budget_closes(ws)
        times, aq = read_array(ws, "gwt.strandaq.bin", "STRAND-AQUEOUS")
        _, srb = read_array(ws, "gwt.strandsrb.bin", "STRAND-SORBED")
        assert np.allclose(aq, mass0 * np.exp(-lam_aq * times), rtol=1e-10), aq
        assert np.allclose(srb, msrb0 * np.exp(-lam_srb * times), rtol=1e-10), srb

    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


def test_fmi(function_tmpdir, targets):
    # the saturation the simulation starts from cannot be read from a budget
    # file, so the water table is held still for the first stress period
    quiet = [(1.0, 0.0)]
    tdis = [(1.0, 1, 1.0), (refill[0][0], 9, 1.0)]

    def build(test):
        ws = test.workspace
        # the flow model runs first, in its own simulation
        fsim = flopy.mf6.MFSimulation(
            sim_name="flow", sim_ws=ws / "flow", exe_name=targets["mf6"]
        )
        flopy.mf6.ModflowTdis(fsim, time_units="DAYS", nper=2, perioddata=tdis)
        add_flow(fsim, strt0, quiet + refill)
        fsim.write_simulation(silent=True)
        success, _ = fsim.run_simulation(silent=True)
        assert success, "flow model failed"

        sim = flopy.mf6.MFSimulation(sim_name="cell", sim_ws=ws, exe_name="mf6")
        flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=2, perioddata=tdis)

        def add_ims(fname):
            return flopy.mf6.ModflowIms(
                sim,
                complexity="moderate",
                outer_dvclose=1e-10,
                inner_dvclose=1e-12,
                linear_acceleration="bicgstab",
                filename=fname,
            )

        add_transport(
            sim,
            add_ims,
            0.0,
            {"stranded_aqueous": mass0},
            fmi=[
                ("GWFHEAD", os.path.join("flow", "gwf.hds")),
                ("GWFBUDGET", os.path.join("flow", "gwf.cbc")),
            ],
        )
        return sim

    def check(test):
        assert_budget_closes(test.workspace)
        check_refill(test.workspace, t0=quiet[0][0])

    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


def test_restart(function_tmpdir, targets):
    half = [(0.9, qin)]
    ws = function_tmpdir

    def run(path, periods, mst_kwargs, strt=strt0, cinit=0.0):
        sim = get_model(path, periods, 9, mst_kwargs, strt=strt, cinit=cinit)
        sim.exe_name = targets["mf6"]
        sim.write_simulation(silent=True)
        success, _ = sim.run_simulation(silent=True)
        assert success, f"simulation in {path} failed"

    # one simulation over both halves, and the first half alone
    run(ws / "full", [half[0], half[0]], {"stranded_aqueous": mass0})
    run(ws / "first", half, {"stranded_aqueous": mass0})

    # the second half starts from what the first half wrote
    hf = flopy.utils.HeadFile(ws / "first" / "gwf.hds")
    head = hf.get_data(totim=hf.get_times()[-1]).flatten()[0]
    _, conc = read_array(ws / "first", "gwt.ucn", "CONCENTRATION")
    _, held = read_array(ws / "first", "gwt.strandaq.bin", "STRAND-AQUEOUS")
    run(
        ws / "second",
        half,
        {"stranded_aqueous": held[-1]},
        strt=head,
        cinit=conc[-1],
    )

    _, cfull = read_array(ws / "full", "gwt.ucn", "CONCENTRATION")
    _, hfull = read_array(ws / "full", "gwt.strandaq.bin", "STRAND-AQUEOUS")
    _, csecond = read_array(ws / "second", "gwt.ucn", "CONCENTRATION")
    _, hsecond = read_array(ws / "second", "gwt.strandaq.bin", "STRAND-AQUEOUS")
    assert np.isclose(csecond[-1], cfull[-1], rtol=1e-6), (
        f"concentration {csecond[-1]} expected {cfull[-1]}"
    )
    assert np.isclose(hsecond[-1], hfull[-1], rtol=1e-6), (
        f"stranded mass {hsecond[-1]} expected {hfull[-1]}"
    )


error_cases = {
    "nooption": (
        {
            "stranded_mass": False,
            "stranded_filerecord": None,
            "stranded_aq_filerecord": None,
            "stranded_srb_filerecord": None,
            "stranded_aqueous": mass0,
        },
        strt0,
        "STRANDED_MASS keyword was not specified",
    ),
    "nosorption": (
        {"stranded_sorbed": mass0},
        strt0,
        "STRANDED_SORBED was specified but sorption is not active",
    ),
    "negative": (
        {"stranded_aqueous": -mass0},
        strt0,
        "cannot be negative",
    ),
    "saturated": (
        {"stranded_aqueous": mass0},
        top + 1.0,
        "is fully saturated at the start of the simulation",
    ),
}


@pytest.mark.parametrize("case", list(error_cases))
def test_errors(function_tmpdir, targets, case):
    mst_kwargs, strt, message = error_cases[case]

    def build(test):
        return get_model(test.workspace, refill, 9, mst_kwargs, strt=strt)

    def check(test):
        listing = (test.workspace / "mfsim.lst").read_text()
        assert message in " ".join(listing.split()), f"expected {message!r}"

    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        xfail=True,
    ).run()

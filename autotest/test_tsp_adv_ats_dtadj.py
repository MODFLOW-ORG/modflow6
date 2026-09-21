"""
Test that a stable time step length submitted to the ATS Package is honored
when DTADJ is zero or one.

DTADJ is the multiplier that ATS applies to the time step in response to the
solver iteration count, and the documentation states that a value of zero or
one means it has no effect on the simulation.  A step length submitted
directly by a package for its own stability (the advective Courant limit from
the ADV Package, or the kinematic wave limit from the SFR Package) is not
scaled by DTADJ and must be honored whatever DTADJ is set to.  Previously the
whole submission was discarded when DTADJ was zero or one, so ATS_PERCEL was
accepted and then silently ignored.

The problem is the 1D column of test_tsp_adv_ats: a well injecting 1.0 into
the first of 100 unit cells with porosity 0.1, so the water crosses one cell
in 0.1 days and ATS_PERCEL = 0.5 limits the step to 0.05 days.  Each case
sets DTADJ to zero or one and starts the period with a DT0 that is either
below the limit, so the step has to grow to it, or above it, so the step has
to shrink to it.  ATS only applies a submitted step from the second step of a
period onward, so the first step is always DT0.  The exact number of steps
follows from that, and it differs from the number taken when the submission
is ignored (500 steps of 0.01 days, or 5 steps of 1.0 days).
"""

import math
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

# name: (transport model type, dtadj, dt0)
cases = {
    "advdtadj0gwt": ("gwt", 0.0, 0.01),
    "advdtadj1gwt": ("gwt", 1.0, 1.0),
    "advdtadj0gwe": ("gwe", 0.0, 1.0),
    "advdtadj1gwe": ("gwe", 1.0, 0.01),
}

nlay, nrow, ncol = 1, 1, 100
perlen = 5.0
delr = 1.0
delc = 1.0
top = 1.0
botm = [0.0]
strt = 1.0
hk = 1.0
porosity = 0.1
qwell = 1.0
percel = 0.5
dtmin = 1.0e-5
rhos = 2700.0
rhow = 1000.0
Cps = 703.7
Cpw = 4183.0

# stable step length submitted by the ADV Package
dtlimit = percel * delr * delc * (top - botm[0]) * porosity / qwell

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def expected_nstp(dt0):
    # first step is dt0, then steps of dtlimit until the end of the period
    return 1 + math.ceil((perlen - dt0) / dtlimit - 1e-9)


def ims(sim, name, acceleration):
    return flopy.mf6.ModflowIms(
        sim,
        ats_outer_maximum_fraction=0.0,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=acceleration,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{name}.ims",
    )


def build_models(name, test):
    mtype, dtadj, dt0 = cases[name]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, 1, 1.0)]
    )
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, dt0, dtmin, perlen, dtadj, 5.0)],
        filename=f"{name}.ats",
    )

    # flow model: steady uniform flow along the column
    gwfname = f"gwf_{name}"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    sim.register_ims_package(ims(sim, gwfname, "CG"), [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, save_specific_discharge=True)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, ncol - 1), 0.0]]}, pname="CHD-1"
    )
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={0: [[(0, 0, 0), qwell, 1.0]]},
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # transport model: upstream weighting with ATS_PERCEL
    tname = f"{mtype}_{name}"
    if mtype == "gwt":
        tsp = flopy.mf6.ModflowGwt(sim, modelname=tname, save_flows=True)
        sim.register_ims_package(ims(sim, tname, "BICGSTAB"), [tsp.name])
        flopy.mf6.ModflowGwtdis(
            tsp,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGwtic(tsp, strt=0.0)
        flopy.mf6.ModflowGwtadv(tsp, scheme="upstream", ats_percel=percel)
        flopy.mf6.ModflowGwtmst(tsp, porosity=porosity)
        flopy.mf6.ModflowGwtssm(tsp, sources=[("WEL-1", "AUX", "CONCENTRATION")])
        flopy.mf6.ModflowGwtoc(
            tsp,
            concentration_filerecord=f"{tname}.ucn",
            saverecord=[("CONCENTRATION", "ALL")],
        )
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=gwfname,
            exgmnameb=tname,
            filename=f"{name}.gwfgwt",
        )
    else:
        tsp = flopy.mf6.ModflowGwe(sim, modelname=tname, save_flows=True)
        sim.register_ims_package(ims(sim, tname, "BICGSTAB"), [tsp.name])
        flopy.mf6.ModflowGwedis(
            tsp,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGweic(tsp, strt=0.0)
        flopy.mf6.ModflowGweadv(tsp, scheme="upstream", ats_percel=percel)
        flopy.mf6.ModflowGweest(
            tsp,
            porosity=porosity,
            heat_capacity_water=Cpw,
            density_water=rhow,
            heat_capacity_solid=Cps,
            density_solid=rhos,
        )
        flopy.mf6.ModflowGwessm(tsp, sources=[("WEL-1", "AUX", "CONCENTRATION")])
        flopy.mf6.ModflowGweoc(
            tsp,
            temperature_filerecord=f"{tname}.ucn",
            saverecord=[("TEMPERATURE", "ALL")],
        )
        flopy.mf6.ModflowGwfgwe(
            sim,
            exgtype="GWF6-GWE6",
            exgmnamea=gwfname,
            exgmnameb=tname,
            filename=f"{name}.gwfgwe",
        )

    return sim, None


def check_output(name, test):
    mtype, dtadj, dt0 = cases[name]
    tname = f"{mtype}_{name}"
    text = "CONCENTRATION" if mtype == "gwt" else "TEMPERATURE"

    fpth = os.path.join(test.workspace, f"{tname}.ucn")
    times = np.array(flopy.utils.HeadFile(fpth, precision="double", text=text).times)
    dt = np.diff(np.concatenate(([0.0], times)))

    # ATS applies the first step as given by DT0
    assert np.isclose(dt[0], dt0), f"first step {dt[0]} is not DT0 {dt0}"

    # every later step must be limited to the value submitted by ADV
    assert np.all(dt[1:] <= dtlimit * (1.0 + 1e-8)), (
        f"steps exceed the ADV limit {dtlimit}: max {dt[1:].max()} with DTADJ {dtadj}"
    )

    # and the step should be raised to the limit as well as lowered to it
    nstp = expected_nstp(dt0)
    assert len(times) == nstp, (
        f"took {len(times)} time steps, expected {nstp} with DTADJ {dtadj} "
        f"and DT0 {dt0}"
    )


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(name, t),
        check=lambda t: check_output(name, t),
    )
    test.run()

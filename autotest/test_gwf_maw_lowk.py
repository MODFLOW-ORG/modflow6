"""
Tests for a multi-aquifer well (MAW) with a single connection to a low
hydraulic conductivity cell, solved with the Newton-Raphson formulation. The
well head is very sensitive in this case, so without help the Newton head can
oscillate and never settle; the MAW package damps these oscillations in
maw_nur.

Two cases are tested:

  1. "maw_lowk_conv" - a feasible well. RATE_SCALING reduces the rate as the
     head drops, so a steady answer exists. The model should converge with no
     over-demand warning.

  2. "maw_lowk_fail" - an infeasible well. A fixed rate asks for more water
     than the cell can supply with no rate reduction, so no steady answer
     exists. The model should fail to converge and print a warning naming the
     well.
"""

import flopy
import pytest
from framework import TestFramework

cases = ["maw_lowk_conv", "maw_lowk_fail"]

# whether each case is expected to fail to converge
xfails = [False, True]

# shared model setup
nlay, nrow, ncol = 1, 1, 11
delr = delc = 100.0
top = 0.0
botm = -50.0
hk = 1.0e-4  # low aquifer hydraulic conductivity
well_col = 5  # interior cell that hosts the single-connection well
well_bottom = -50.0
pump_rate = -20.0  # extraction (negative), larger than the cell can supply
pump_elev = -45.0
reduction_length = 5.0


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

    # Use no solver-level under-relaxation so the test exercises the MAW
    # package's own head damping (maw_nur), which is turned on by the model
    # NEWTON UNDER_RELAXATION option.
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=1e-7,
        outer_maximum=200,
        under_relaxation="NONE",
        inner_maximum=100,
        inner_dvclose=1e-8,
        linear_acceleration="BICGSTAB",
    )

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON UNDER_RELAXATION"
    )
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
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=hk)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, steady_state={0: True})

    # constant heads at both ends supply water to the aquifer
    chd = [[(0, 0, 0), 0.0], [(0, 0, ncol - 1), 0.0]]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)

    # single-connection MAW well
    perioddata = [(0, "rate", pump_rate)]
    if not xfails[idx]:
        # feasible case: reduce the rate as the head drops so an answer exists
        perioddata.append((0, "rate_scaling", pump_elev, reduction_length))
    flopy.mf6.ModflowGwfmaw(
        gwf,
        print_head=True,
        no_well_storage=True,
        packagedata=[[0, 0.15, well_bottom, 0.0, "THIEM", 1]],
        connectiondata=[[0, 0, (0, 0, well_col), 0.0, well_bottom, 0.0, 0.0]],
        perioddata={0: perioddata},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim


def check_output(idx, test):
    with open(test.workspace / "mfsim.lst") as f:
        listing = f.read()

    warning = "aquifer can supply" in listing

    if xfails[idx]:
        # the infeasible well should trigger the over-demand warning, and the
        # warning should name the well
        assert warning, "expected the MAW over-demand warning, but it was not found"
        assert "MAW well" in listing
        assert "(MAW_1-1)" in listing, "warning did not identify the well"
    else:
        # the feasible well should converge with no over-demand warning
        assert "Normal termination" in listing, "model did not terminate normally"
        assert not warning, "unexpected MAW over-demand warning for a feasible well"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        xfail=xfails[idx],
    )
    test.run()

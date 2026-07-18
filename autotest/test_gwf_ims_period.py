"""
Tests the stress-period-varying IMS linear settings (LINEAR_PERIODDATA).

The feature is not yet supported directly by flopy, so the base model is built
with flopy and the LINEAR_PERIODDATA option plus the period-data file are
injected into the written input.

Cases:
  - ims_period       : PERIOD 2 raises INNER_MAXIMUM (which grows the solver
                       convergence-history arrays, exercising their runtime
                       reallocation) and switches the preconditioner ILU0 -> ILUT
                       (forces the preconditioner to be reallocated); PERIOD 3
                       lowers INNER_MAXIMUM again and tightens INNER_DVCLOSE. Must
                       run and process both period blocks.
  - ims_period_badkw : a keyword not allowed in a PERIOD block (SCALING_METHOD)
                       must terminate with an error (xfail).
"""

import os
import re

import flopy
import pytest
from framework import TestFramework

cases = ["ims_period", "ims_period_badkw"]
xfail = [False, True]

nlay, nrow, ncol = 1, 5, 5
nper = 3


def _inject_period_data(ws, name, period_text):
    # flopy writes "BEGIN options" (mixed case), so match case-insensitively
    ims = os.path.join(ws, f"{name}.ims")
    with open(ims) as f:
        txt = f.read()
    new = re.sub(
        r"(?i)(BEGIN\s+OPTIONS[^\n]*\n)",
        rf"\1  LINEAR_PERIODDATA FILEIN {name}.impd\n",
        txt,
        count=1,
    )
    assert new != txt, "could not inject LINEAR_PERIODDATA into the IMS OPTIONS block"
    with open(ims, "w") as f:
        f.write(new)
    with open(os.path.join(ws, f"{name}.impd"), "w") as f:
        f.write(period_text)


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim, complexity="simple", print_option="all")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=10.0, delc=10.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 10.0], [(0, nrow - 1, ncol - 1), 5.0]],
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation()

    if idx == 0:
        period = (
            "BEGIN PERIOD 2\n"
            "  INNER_MAXIMUM 500\n"
            "  PRECONDITIONER_LEVELS 5\n"
            "  PRECONDITIONER_DROP_TOLERANCE 1.0e-4\n"
            "END PERIOD\n\n"
            "BEGIN PERIOD 3\n"
            "  INNER_MAXIMUM 50\n"
            "  INNER_DVCLOSE 1.0e-9\n"
            "END PERIOD\n"
        )
    else:
        period = "BEGIN PERIOD 2\n  SCALING_METHOD L2NORM\nEND PERIOD\n"
    _inject_period_data(ws, name, period)
    return sim, None


def check_output(idx, test):
    lst = os.path.join(test.workspace, "mfsim.lst")
    with open(lst) as f:
        text = f.read()
    assert "LINEAR_PERIODDATA input will be read" in text, (
        "the period-data file was not opened"
    )
    for kper in (2, 3):
        assert f"PROCESSING LINEAR PERIOD DATA FOR STRESS PERIOD {kper}" in text, (
            f"linear settings for stress period {kper} were not processed"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=(lambda t: check_output(idx, t)) if idx == 0 else None,
        xfail=xfail[idx],
        overwrite=False,  # keep the injected LINEAR_PERIODDATA + period file
    )
    test.run()

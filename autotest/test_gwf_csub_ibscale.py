"""
Tests scaling of interbed thicknesses that exceed the cell thickness.

Interbed thicknesses in a cell sum to slightly more than the cell thickness, so
the summed interbed thickness exceeds it by a numerical-roundoff amount. CSUB
proportionally scales the cell's interbeds to fit instead of terminating. Each
cell has two interbeds (one cell with two no-delay interbeds, one with a
no-delay and a delay interbed) so the summed excess is exercised.

Cases:
  - csub_ibscale    : interbeds specified as cell fractions; both cells scaled.
  - csub_ibscale_th : same, but interbeds specified as absolute thicknesses.
  - csub_ibscale_er : the summed interbed thickness exceeds the cell thickness
                      by a gross amount, which is an error (xfail).
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_ibscale", "csub_ibscale_th", "csub_ibscale_er"]
use_frac = [True, False, True]  # csub_ibscale_th uses absolute thicknesses
gross = [False, False, True]

# cell thickness is top - botm = 10; the roundoff tolerance is 1e-5 * 10 = 1e-4
top, botm = 10.0, 0.0
# fractions that sum to just over one (overshoot 5e-5 < 1e-4 tolerance)
frac_a, frac_b = 0.5, 0.500005
# a gross overshoot (0.5 well beyond the tolerance)
frac_gross = 0.6


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    baq = top - botm
    fb = frac_gross if gross[idx] else frac_b
    # interbed thickness column: cell fractions or absolute thicknesses
    if use_frac[idx]:
        col_a, col_b = frac_a, fb
    else:
        col_a, col_b = frac_a * baq, fb * baq
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0e3, 10, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple", outer_dvclose=1e-8)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    # cell 0 drives the two active csub cells 1 and 2
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=top)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    # CSUB provides all storage, so STO specific storage must be zero
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=0, ss=0.0, sy=0.0, steady_state={0: True}, transient={1: True}
    )
    chd = flopy.mf6.ModflowGwfchd(
        gwf, maxbound=1, stress_period_data={0: [[(0, 0, 0), "driver"]]}
    )
    chd.ts.initialize(
        filename=f"{name}.ts",
        timeseries=[(0.0, top), (1.0, top), (1.0 + 1.0e3, top - 5.0)],
        time_series_namerecord=["driver"],
        interpolation_methodrecord=["linear"],
    )
    # icsubno, cellid, cdelay, pcs0, thick, rnb, ssv_cc, sse_cr, theta, kv, h0
    pkgdata = [
        # cell 1: two no-delay interbeds summing to just over the cell thickness
        (0, (0, 0, 1), "nodelay", 0.0, col_a, 1.0, 1e-4, 1e-5, 0.2, 1e-3, 0.0),
        (1, (0, 0, 1), "nodelay", 0.0, col_b, 1.0, 1e-4, 1e-5, 0.2, 1e-3, 0.0),
        # cell 2: a no-delay and a delay interbed summing to just over it
        (2, (0, 0, 2), "nodelay", 0.0, col_a, 1.0, 1e-4, 1e-5, 0.2, 1e-3, 0.0),
        (3, (0, 0, 2), "delay", 0.0, col_b, 1.0, 1e-4, 1e-5, 0.2, 1e-3, 0.0),
    ]
    flopy.mf6.ModflowGwfcsub(
        gwf,
        cell_fraction=use_frac[idx],
        ninterbeds=4,
        ndelaycells=7,
        gammaw=9806.65,
        beta=0.0,
        sgm=1.7,
        sgs=2.0,
        cg_theta=0.2,
        cg_ske_cr=1e-5,
        packagedata=pkgdata,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation()
    return sim, None


def check_output(idx, test):
    # only the roundoff (scaled) case is checked; the gross case is xfail.
    # the scaling is silent, so it is verified by the run completing normally
    # with finite heads and a closed budget
    name = cases[idx]
    with open(os.path.join(test.workspace, f"{name}.lst")) as f:
        text = f.read()

    heads = flopy.utils.HeadFile(
        os.path.join(test.workspace, f"{name}.hds")
    ).get_alldata()
    assert np.isfinite(heads).all(), "non-finite heads"
    cum = [
        abs(float(m.group(1)))
        for m in re.finditer(
            r"PERCENT DISCREPANCY =\s*([-0-9.E+]+)\s+PERCENT DISCREPANCY", text
        )
    ]
    assert cum, "no cumulative budget discrepancy reported"
    assert max(cum) < 0.1, f"mass balance did not close (max discrepancy {max(cum)}%)"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=(lambda t: check_output(idx, t)) if not gross[idx] else None,
        xfail=gross[idx],
    )
    test.run()

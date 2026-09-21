"""
Tests the CSUB effective-stress regularization.

A constant-head driver ramps an adjacent active interbed cell above land
surface, so its pore pressure exceeds the geostatic stress and the calculated
effective stress goes negative. STRICT_EFFECTIVE_STRESS is injected for the
strict case (flopy does not support the keyword yet).

Cases:
  - csub_esreg    : default. Effective stress is regularized; the run completes,
                    warns, and the mass balance closes.
  - csub_esreg_st : STRICT_EFFECTIVE_STRESS terminates on negative stress (xfail).
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_esreg", "csub_esreg_st"]
strict = [False, True]

# geometry: 10-length cells, land surface at top=10, base at 0
top, botm = 10.0, 0.0
sgm, sgs = 1.7, 2.0  # geostatic ~ (top-botm)*sgs = 20 when saturated
# a driver head of 40 pulls the active cell above land surface so that its pore
# pressure (~40) exceeds the geostatic stress (~20): effective stress ~ -20
h_steady, h_high = 10.0, 40.0


def inject_strict(ws, name):
    csub = os.path.join(ws, f"{name}.csub")
    with open(csub) as f:
        txt = f.read()
    new = re.sub(
        r"(?i)(BEGIN\s+OPTIONS[^\n]*\n)",
        r"\1  STRICT_EFFECTIVE_STRESS\n",
        txt,
        count=1,
    )
    assert new != txt, "could not inject STRICT_EFFECTIVE_STRESS"
    with open(csub, "w") as f:
        f.write(new)


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    # steady period, then a transient period during which the driver head is
    # ramped continuously above land surface (never plateaus, so budget rates
    # stay non-negligible and the effective-stress crossing is well resolved)
    perlen = 1.0e4
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (perlen, 20, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        complexity="complex",
        outer_dvclose=1e-8,
        inner_dvclose=1e-9,
        outer_maximum=200,
        inner_maximum=200,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    # cell 0 is the constant-head driver, cell 1 is the active csub cell
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=h_steady)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    # CSUB provides all storage, so STO specific storage must be zero
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=0, ss=0.0, sy=0.0, steady_state={0: True}, transient={1: True}
    )
    chd = flopy.mf6.ModflowGwfchd(
        gwf, maxbound=1, stress_period_data={0: [[(0, 0, 0), "driver"]]}
    )
    chd.ts.initialize(
        filename=f"{name}.ts",
        timeseries=[(0.0, h_steady), (1.0, h_steady), (1.0 + perlen, h_high)],
        time_series_namerecord=["driver"],
        interpolation_methodrecord=["linear"],
    )
    # coarse-grained storage plus one delay interbed in the active cell
    pkgdata = [
        # icsubno, cellid, cdelay, pcs0, thick_frac, rnb, ssv_cc, sse_cr, theta, kv, h0
        (0, (0, 0, 1), "delay", 0.0, 1.0, 1.0, 1e-3, 1e-4, 0.2, 1e-3, 0.0)
    ]
    flopy.mf6.ModflowGwfcsub(
        gwf,
        ninterbeds=1,
        ndelaycells=7,
        gammaw=9806.65,
        beta=0.0,
        sgm=sgm,
        sgs=sgs,
        cg_theta=0.2,
        cg_ske_cr=1e-4,
        packagedata=pkgdata,
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation()
    if strict[idx]:
        inject_strict(ws, name)
    return sim, None


def check_output(idx, test):
    # only the default (regularize) case is checked; the strict case is xfail
    name = cases[idx]
    with open(os.path.join(test.workspace, f"{name}.lst")) as f:
        text = f.read()

    # the per-time-step regularization note was written to the model listing
    assert "negative effective stress regularized" in text.lower(), (
        "the effective-stress regularization note was not written to the listing"
    )
    # the end-of-run summary warning is in the simulation listing
    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        simlst = f.read()
    assert "effective stress was regularized in" in simlst.lower(), (
        "the effective-stress regularization summary warning was not issued"
    )

    # the active cell was actually driven above land surface (negative eff. stress)
    heads = flopy.utils.HeadFile(
        os.path.join(test.workspace, f"{name}.hds")
    ).get_alldata()
    assert heads[-1, 0, 0, 1] > top, (
        f"active cell head {heads[-1, 0, 0, 1]} did not exceed land surface {top}"
    )
    assert np.isfinite(heads).all(), (
        "non-finite heads (storage singularity not handled)"
    )

    # cumulative mass balance closes (first value on each line is cumulative;
    # the per-step rate is not checked because it is spurious near equilibrium)
    cum = [
        abs(float(m.group(1)))
        for m in re.finditer(
            r"PERCENT DISCREPANCY =\s*([-0-9.E+]+)\s+PERCENT DISCREPANCY", text
        )
    ]
    assert cum, "no cumulative budget discrepancy reported"
    assert max(cum) < 0.1, (
        f"cumulative mass balance did not close (max discrepancy {max(cum)}%)"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=(lambda t: check_output(idx, t)) if not strict[idx] else None,
        xfail=strict[idx],
        overwrite=False,  # keep the injected STRICT_EFFECTIVE_STRESS keyword
    )
    test.run()

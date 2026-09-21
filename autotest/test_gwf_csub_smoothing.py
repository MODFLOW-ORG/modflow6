"""
Feature test for the CSUB ELASTIC_INELASTIC_SMOOTHING option.

A small confined delay-interbed problem is drawn down past the preconsolidation
stress and run twice, with the option off and on. With the Newton formulation
the delay interbed skeletal storage switches abruptly from elastic to inelastic
at the preconsolidation stress; the option smooths that switch over a small
effective-stress window.

Tests that the option is read and shifts delay interbed storage from inelastic
toward elastic near the switch: with smoothing the elastic delay interbed
storage goes from zero to nonzero and the inelastic storage decreases. The
reported elastic/inelastic split is weighted by the smoothing fraction, so it
only repartitions storage between the two budget terms and the total delay
interbed storage is conserved. Also confirms the smoothed storage is applied
consistently in the budget by checking that the mass balance closes (the
smoothed ssk is used in both the solve and the storage budget).

Does not test the convergence improvement the option provides. That failure is
an emergent property of a large coupled model and cannot be reproduced in a
small fast model; it is verified on the model that motivated the option (see
MODFLOW-USGS/modflow6 issue #2896).
"""

import os
import re

import flopy
import pytest
from framework import TestFramework

cases = ["csub_smoothing"]
smooth_ws = "smooth"  # subdirectory for the smoothing-on simulation

# spatial discretization: layer 1 holds the delay interbed, layer 2 is the
# specified-head driver; deep so the preconsolidation stress (and the smoothing
# window) are realistic
nlay, nrow, ncol = 2, 1, 1
top = 0.0
botm = [-500.0, -1000.0]

# temporal discretization: steady first period, then a 30-year drawdown
nper = 2
perlen = [1.0, 10957.5]
nstp = [1, 40]
hdraw = -120.0  # layer 2 head at the end of the transient period

# properties (1a Chicot interbed: Sskv=2e-5, Sskv:Sske=100, normally consolidated)
ssv, sse, theta = 2e-5, 2e-7, 0.3
dbthick, kv = 200.0, 1e-4
gammaw, beta, sgm, sgs = 62.4, 2.227e-8, 1.7, 2.0


def build_sim(ws, name, smoothing):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    tdis_rc = [(perlen[i], nstp[i], 1.0) for i in range(nper)]
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-6,
        outer_maximum=100,
        inner_dvclose=1e-6,
        inner_maximum=100,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="NEWTON")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1000.0,
        delc=1000.0,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0, k33=2e-3)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=0,
        ss=0.0,
        sy=0.0,
        steady_state={0: True},
        transient={1: True},
    )
    # layer 2 specified head, ramped down over the transient period
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=1,
        stress_period_data={0: [[(1, 0, 0), "ct"]]},
    )
    chd.ts.initialize(
        filename=f"{name}.ts",
        timeseries=[(0.0, 0.0), (perlen[0], 0.0), (perlen[0] + perlen[1], hdraw)],
        time_series_namerecord=["ct"],
        interpolation_methodrecord=["linear"],
    )
    # a single delay interbed in the free layer (layer 1)
    pkgdata = [(0, (0, 0, 0), "delay", 0.0, dbthick, 1.0, ssv, sse, theta, kv, 0.0)]
    flopy.mf6.ModflowGwfcsub(
        gwf,
        ndelaycells=19,
        ninterbeds=1,
        gammaw=gammaw,
        beta=beta,
        sgm=sgm,
        sgs=sgs,
        cg_theta=theta,
        cg_ske_cr=sse,
        elastic_inelastic_smoothing=smoothing,
        packagedata=pkgdata,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    # same model with the smoothing option off (base) and on (comparison)
    sim = build_sim(test.workspace, name, smoothing=False)
    mc = build_sim(os.path.join(test.workspace, smooth_ws), name, smoothing=True)
    return sim, mc


def csub_storage(ws, name):
    # cumulative elastic and inelastic delay interbed storage from the listing
    fpth = os.path.join(ws, f"{name}.lst")
    bud = flopy.utils.Mf6ListBudget(fpth).get_cumulative()
    inelastic = float(bud["CSUB-INELASTIC_IN"][-1])
    elastic = float(bud["CSUB-ELASTIC_IN"][-1])
    return inelastic, elastic


def max_discrepancy(ws, name):
    # largest absolute percent discrepancy reported in the listing
    fpth = os.path.join(ws, f"{name}.lst")
    with open(fpth) as f:
        vals = re.findall(r"PERCENT DISCREPANCY =\s*([-0-9.E+]+)", f.read())
    return max(abs(float(v)) for v in vals)


def check_output(idx, test):
    name = cases[idx]
    off_ine, off_ela = csub_storage(test.workspace, name)
    on_ine, on_ela = csub_storage(os.path.join(test.workspace, smooth_ws), name)

    # the drawdown drives the interbed inelastic in both runs (switch is active)
    assert off_ine > 1.0e3, f"interbed did not load inelastically ({off_ine})"
    # without smoothing the hard switch leaves no near-switch elastic storage
    assert abs(off_ela) < 1.0, (
        f"unexpected elastic storage without smoothing ({off_ela})"
    )
    # with smoothing, storage near the switch shifts from inelastic to elastic
    assert on_ela > 1.0, f"smoothing did not add elastic storage ({on_ela})"
    assert on_ine < off_ine, (
        f"smoothing did not reduce inelastic storage ({on_ine} !< {off_ine})"
    )
    # the weighted split only repartitions storage between the elastic and
    # inelastic budget terms, so the total delay interbed storage is conserved
    # aside from the small in-window blend of the skeletal storage coefficient
    off_tot, on_tot = off_ine + off_ela, on_ine + on_ela
    assert abs(on_tot - off_tot) < 0.01 * off_tot, (
        f"smoothing changed total delay interbed storage beyond the blend "
        f"({on_tot} vs {off_tot})"
    )

    # smoothed ssk must be applied consistently in the budget -> mass balance closes
    for ws in (test.workspace, os.path.join(test.workspace, smooth_ws)):
        pd = max_discrepancy(ws, name)
        assert pd < 1.0e-3, f"nonzero mass-balance discrepancy ({pd}) in {ws}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

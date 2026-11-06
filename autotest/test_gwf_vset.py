"""
Tests for the variably-saturated evapotranspiration package VSET.

Case "vset-feddes": a 1D unsaturated column (VSF/Richards) is driven to steady
state with a water table fixed near mid-column by a constant-head cell at the
base.  The VSET package applies a Feddes pressure-head root-water uptake to every
active cell.  Because the simulated pressure head varies with depth, the cells
span the Feddes stress regimes (anoxia, wet limb, optimal, dry limb, wilting).
The check verifies that each cell's VSET budget flow equals the Feddes reduction
applied to the *simulated* pressure head, i.e.

    q_vset = -alpha(psi) * rate * area ,   psi = h - z_cellcenter

which is an exact self-consistency check of the uptake kernel, the budget
output, area weighting, and sign convention.

Case "vset-no-vsf": the same VSET package without an active VSF package.  VSET may
only be applied on VSF (Richards) cells, so the run must terminate with an
error.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from gwf_test_utils import get_vsf_soil_data

cases = ["vset-feddes", "vset-no-vsf"]
with_vsf = [True, False]
xfail = [False, True]

# Feddes pressure-head thresholds (cm), h1 > h2 > h3 > h4
H1 = -2.0
H2 = -8.0
H3 = -25.0
H4 = -40.0

# per-cell potential transpiration rate (cm/s); small to keep the column
# close to hydrostatic so the regimes are well separated with depth
POT_RATE = 1.0e-7

DELR = 10.0
DELC = 10.0


def feddes_alpha(psi):
    if psi >= H1 or psi <= H4:
        return 0.0
    elif psi >= H2:
        return (psi - H1) / (H2 - H1)
    elif psi >= H3:
        return 1.0
    else:
        return (psi - H4) / (H3 - H4)


def build_models(idx, test):
    nlay, nrow, ncol = 20, 1, 1
    column_height = 100.0
    delz = column_height / nlay
    top = column_height
    botm = [column_height - (k + 1) * delz for k in range(nlay)]
    hk = 0.00944  # cm/s

    # water table near mid-column via constant head at the base
    h_base = 50.0
    strt = np.full((nlay, nrow, ncol), h_base, dtype=float)

    nouter, ninner = 200, 300
    hclose, rclose, relax = 1e-7, 1e-7, 0.97

    name = cases[idx]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="SECONDS", nper=1, perioddata=[(1.0, 1, 1.0)])

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam", save_flows=True
    )

    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.9,
        under_relaxation_kappa=0.0001,
        under_relaxation_gamma=0.0,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=DELR,
        delc=DELC,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")
    flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, icelltype=0, k=hk, k33=hk
    )
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-4, sy=0.2, iconvert=0, steady_state=True)

    # unsaturated zone Richards flow (omitted for the negative case)
    if with_vsf[idx]:
        soil = get_vsf_soil_data("Celia1990-eq10-Haverkamp")
        flopy.mf6.ModflowGwfvsf(
            gwf,
            iunsat=1,
            storage_scheme="chord-slope",
            kr_averaging="harmonic",
            soil_model="Haverkamp",
            porosity=soil["porosity"],
            satres=soil["satres"],
            alphahvk=soil["alpha"],
            nhvk=soil["n"],
            betahvk=soil["beta"],
            khvk=soil["k"],
        )

    # variably-saturated evapotranspiration (Feddes), applied to all cells
    # except the constant-head base cell
    spd = [[(k, 0, 0), POT_RATE] for k in range(nlay - 1)]
    flopy.mf6.ModflowGwfvset(
        gwf,
        et_model="feddes",
        h1=H1,
        h2=H2,
        h3=H3,
        h4=H4,
        maxbound=len(spd),
        stress_period_data={0: spd},
        print_flows=True,
        save_flows=True,
        pname="VSET-1",
    )

    # constant head at base sets the water table
    c = {0: [[(nlay - 1, 0, 0), h_base]]}
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=c, pname="CHD-1")

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_feddes(idx, test):
    model_name = "gwf_" + cases[idx]
    ws = test.workspace

    # simulated heads
    hds = flopy.utils.HeadFile(os.path.join(ws, f"{model_name}.hds"))
    head = hds.get_data()[:, 0, 0]

    # cell-center elevations
    grb = flopy.mf6.utils.MfGrdFile(os.path.join(ws, f"{model_name}.dis.grb"))
    mg = grb.modelgrid
    nlay = mg.nlay
    dz = mg.delz.flatten()
    botm = mg.botm.flatten()
    zcenter = np.array([botm[k] + 0.5 * dz[k] for k in range(nlay)])

    area = DELR * DELC

    # VSET cell-by-cell flows
    cbc = flopy.utils.CellBudgetFile(os.path.join(ws, f"{model_name}.cbc"))
    vset = cbc.get_data(text="VSET")[-1]
    q_by_node = {int(rec["node"]): float(rec["q"]) for rec in vset}

    n_partial = 0
    n_optimal = 0
    for k in range(nlay - 1):
        node = k + 1  # 1-based node number for (k, 0, 0) in a 1-col grid
        psi = head[k] - zcenter[k]
        alpha = feddes_alpha(psi)
        expected = -alpha * POT_RATE * area
        assert node in q_by_node, f"missing VSET flow for node {node}"
        sim_q = q_by_node[node]
        assert np.isclose(sim_q, expected, rtol=1e-6, atol=1e-12), (
            f"cell {k}: psi={psi:.3f} alpha={alpha:.4f} "
            f"expected {expected:.6e}, got {sim_q:.6e}"
        )
        if 0.0 < alpha < 1.0:
            n_partial += 1
        elif alpha == 1.0:
            n_optimal += 1

    # make sure the test actually exercises the stress function, not only
    # the trivial alpha == 0 regime
    assert n_optimal >= 1, "no cell reached the optimal (alpha=1) regime"
    assert n_partial >= 1, "no cell reached a partial-stress regime"


def check_no_vsf(idx, test):
    # the run must terminate because VSET requires VSF
    tag = "VSET package requires the VSF package"
    with open(os.path.join(test.workspace, "mfsim.lst"), "r") as f:
        found = any(tag in line for line in f)
    assert found, "expected VSF-required error not found in mfsim.lst"


def check_output(idx, test):
    if with_vsf[idx]:
        check_feddes(idx, test)
    else:
        check_no_vsf(idx, test)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=xfail[idx],
    )
    test.run()

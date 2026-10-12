"""Uniform density must preserve the DRN discharge-depth head solution."""

from itertools import product

import flopy
import numpy as np
import pytest
from framework import TestFramework

settings = list(product([-1.0, 1.0], [False, True], [1.0, 2.0], [0.0, 16.0]))
cases = [f"buydrndepth{i}" for i in range(len(settings))]


def build_models(idx, test):
    depth, newton, multiplier, concentration = settings[idx]
    # Both depth signs describe the same interval [1, 2]. The solution h=1.5
    # has linear factor 1/2 or Newton cubic factor 3/8. Uniform density cancels
    # from the cell balance, so CHD heads below apply to saline and fresh water.
    factor = 0.375 if newton else 0.5
    inflow_head = 1.5 + 0.5 * multiplier * factor
    elevation = 1.0 if depth > 0.0 else 2.0
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="flow",
        save_flows=True,
        newtonoptions="NEWTON" if newton else None,
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        outer_maximum=200,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(ims, ["flow"])
    grid = {
        "nlay": 1,
        "nrow": 1,
        "ncol": 2,
        "delr": 1.0,
        "delc": 1.0,
        "top": 1.0,
        "botm": 0.0,
    }
    flopy.mf6.ModflowGwfdis(gwf, **grid)
    flopy.mf6.ModflowGwfic(gwf, strt=inflow_head + 0.05)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0, icelltype=0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        pname="CHD",
        auxiliary=["CONCENTRATION"],
        stress_period_data=[[(0, 0, 0), inflow_head, concentration]],
    )
    flopy.mf6.ModflowGwfdrn(
        gwf,
        pname="DRN",
        auxiliary=["DEPTH", "MULTIPLIER"],
        auxdepthname="DEPTH",
        auxmultname="MULTIPLIER",
        stress_period_data=[[(0, 0, 1), elevation, 1.0, depth, multiplier]],
    )
    flopy.mf6.ModflowGwfbuy(
        gwf,
        packagedata=[(0, 1.3, 0.0, "trans", "CONCENTRATION")],
        density_filerecord="flow.density",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="flow.hds",
        budget_filerecord="flow.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans", save_flows=True)
    imst = flopy.mf6.ModflowIms(
        sim,
        filename="trans.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-11,
        rcloserecord="1e-9 strict",
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(imst, ["trans"])
    flopy.mf6.ModflowGwtdis(gwt, **grid)
    flopy.mf6.ModflowGwtic(gwt, strt=concentration)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="trans.ucn",
        budget_filerecord="trans.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea="flow", exgmnameb="trans")
    return sim, None


def check_output(idx, test):
    _, newton, multiplier, concentration = settings[idx]
    factor = 0.375 if newton else 0.5
    relative_density = 1.0 + 1.3 * concentration / 1000.0
    inflow_head = 1.5 + 0.5 * multiplier * factor
    expected_discharge = 0.5 * multiplier * factor * relative_density
    ws = test.workspace
    h = flopy.utils.HeadFile(ws / "flow.hds", precision="double").get_data().ravel()
    c = (
        flopy.utils.HeadFile(ws / "trans.ucn", precision="double", text="CONCENTRATION")
        .get_data()
        .ravel()
    )
    rho = (
        flopy.utils.HeadFile(ws / "flow.density", precision="double", text="DENSITY")
        .get_data()
        .ravel()
    )
    np.testing.assert_allclose(h, [inflow_head, 1.5], atol=1e-8, rtol=0.0)
    np.testing.assert_allclose(rho, 1000.0 * relative_density, atol=1e-8, rtol=0.0)
    budget = flopy.utils.CellBudgetFile(ws / "flow.cbc", precision="double")
    qdrn = budget.get_data(text="DRN")[0]["q"].item()
    qchd = budget.get_data(text="CHD")[0]["q"].item()
    qface = budget.get_data(text="FLOW-JA-FACE")[0].ravel()[3]
    np.testing.assert_allclose(
        [qdrn, qchd, qface],
        [-expected_discharge, expected_discharge, expected_discharge],
        atol=1e-8,
        rtol=0.0,
    )
    np.testing.assert_allclose(qdrn + qchd, 0.0, atol=1e-8, rtol=0.0)
    # Exact upwind inventory for the recorded finite-tolerance water rates.
    c0 = concentration * (0.35 + qchd) / (0.35 + qface)
    c1 = (0.35 * concentration + qface * c0) / (0.35 - qdrn)
    np.testing.assert_allclose(c, [c0, c1], atol=1e-8, rtol=0.0)
    mass_residual = [
        0.35 * (c[0] - concentration) - qchd * concentration + qface * c[0],
        0.35 * (c[1] - concentration) - qface * c[0] - qdrn * c[1],
    ]
    np.testing.assert_allclose(mass_residual, 0.0, atol=1e-8, rtol=0.0)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

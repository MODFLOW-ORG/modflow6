"""A weak GHB anchor must not cause a density-dependent RIV iteration cycle."""

import flopy
import numpy as np
import pytest
from framework import TestFramework

settings = [
    (rhs, solver, concentration)
    for rhs, solver in [(True, "CG"), (True, "BICGSTAB"), (False, "BICGSTAB")]
    for concentration in [0.0, 16.0]
]
cases = [f"buyrivrhs{i}" for i in range(len(settings))]


def build_models(idx, test):
    rhs, solver, concentration = settings[idx]
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname="flow", save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        filename="flow.ims",
        outer_maximum=100,
        inner_maximum=100,
        outer_dvclose=1e-10,
        inner_dvclose=1e-10,
        rcloserecord=1e-8,
        linear_acceleration=solver,
    )
    sim.register_ims_package(ims, [gwf.name])
    grid = dict(nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=1.0, botm=0.0)
    flopy.mf6.ModflowGwfdis(gwf, **grid)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfwel(
        gwf,
        pname="WEL",
        auxiliary="CONCENTRATION",
        stress_period_data=[[(0, 0, 0), 0.9, 0.0]],
    )
    flopy.mf6.ModflowGwfghb(
        gwf,
        pname="GHB",
        auxiliary="CONCENTRATION",
        stress_period_data=[[(0, 0, 0), 0.0, 0.01, 0.0]],
    )
    flopy.mf6.ModflowGwfriv(
        gwf,
        pname="RIV",
        auxiliary="CONCENTRATION",
        stress_period_data=[[(0, 0, 0), 1.0, 10.0, 1.0, 0.0]],
    )
    flopy.mf6.ModflowGwfbuy(
        gwf,
        hhformulation_rhs=rhs,
        density_filerecord="flow.density",
        packagedata=[(0, 1.3, 0.0, "transport", "CONCENTRATION")],
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="flow.hds",
        budget_filerecord="flow.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname="transport")
    imst = flopy.mf6.ModflowIms(
        sim,
        filename="transport.ims",
        outer_dvclose=1e-10,
        inner_dvclose=1e-10,
        rcloserecord=1e-8,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(imst, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, **grid)
    flopy.mf6.ModflowGwtic(gwt, strt=concentration)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=[(pname, "AUX", "CONCENTRATION") for pname in ["WEL", "GHB", "RIV"]],
    )
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="transport.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgmnamea=gwf.name, exgmnameb=gwt.name)
    return sim, None


def check_output(idx, test):
    _, _, concentration = settings[idx]
    ws = test.workspace
    rho = 1.0 + 1.3 * concentration / 1000.0
    h = flopy.utils.HeadFile(ws / "flow.hds", precision="double").get_data().item()
    d = (
        flopy.utils.HeadFile(ws / "flow.density", text="DENSITY", precision="double")
        .get_data()
        .item()
    )
    np.testing.assert_allclose(d, rho * 1000.0, rtol=0.0, atol=1e-10)
    # Pressure at z=0.5 is rho*(h-0.5); freshwater boundary pressure is hb-0.5.
    # 0 = Q + Cr*(stage-rho*h+(rho-1)*z) + Cg*(-rho*h+(rho-1)*z).
    expected = (0.9 + 10.0 + 10.01 * (rho - 1.0) * 0.5) / (10.01 * rho)
    assert expected > 1.0, "Analytic river must be on the active branch"
    np.testing.assert_allclose(h, expected, rtol=0.0, atol=1e-8)
    qr = 10.0 * (1.0 - rho * h + (rho - 1.0) * 0.5)
    qg = 0.01 * (-rho * h + (rho - 1.0) * 0.5)
    assert qr < 0.0 and qg < 0.0
    budget = flopy.utils.CellBudgetFile(ws / "flow.cbc", precision="double")
    reported = [
        budget.get_data(text=name)[0]["q"].item() for name in ["WEL", "RIV", "GHB"]
    ]
    np.testing.assert_allclose(reported, [0.9, qr, qg], rtol=0.0, atol=1e-8)
    np.testing.assert_allclose(sum(reported), 0.0, rtol=0.0, atol=1e-8)
    c = (
        flopy.utils.HeadFile(
            ws / "transport.ucn", text="CONCENTRATION", precision="double"
        )
        .get_data()
        .item()
    )
    # Unit volume and time step, porosity .35; inflow has zero concentration.
    np.testing.assert_allclose(
        0.35 * (c - concentration), (qr + qg) * c, rtol=0.0, atol=1e-8
    )


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

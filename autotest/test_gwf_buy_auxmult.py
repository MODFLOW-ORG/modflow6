"""AUXMULTNAME must scale the entire BUY boundary conductance."""

import itertools

import flopy
import numpy as np
import pytest
from framework import TestFramework

settings = list(itertools.product(["GHB", "RIV", "DRN"], [0.0, 2.0], [0.0, 16.0]))
cases = [f"buyauxmult{i}" for i in range(len(settings))]


def build_models(idx, test):
    kind, multiplier, concentration = settings[idx]
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, perioddata=[(1.0, 1, 1.0)])
    grid = dict(nlay=1, nrow=1, ncol=2, top=1.0, botm=0.0, delr=1.0, delc=1.0)
    for suffix, use_multiplier in [("mult", True), ("direct", False)]:
        flow, transport = f"flow_{suffix}", f"gwt_{suffix}"
        gwf = flopy.mf6.ModflowGwf(sim, modelname=flow, save_flows=True)
        ims = flopy.mf6.ModflowIms(
            sim,
            filename=f"{flow}.ims",
            outer_dvclose=1e-10,
            inner_dvclose=1e-10,
            rcloserecord=1e-9,
            linear_acceleration="BICGSTAB",
        )
        sim.register_ims_package(ims, [flow])
        flopy.mf6.ModflowGwfdis(gwf, **grid)
        flopy.mf6.ModflowGwfic(gwf, strt=2.0)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
        flopy.mf6.ModflowGwfchd(
            gwf,
            pname="CHD",
            auxiliary="CONCENTRATION",
            stress_period_data=[[(0, 0, 0), 2.0, concentration]],
        )
        conductance = 10.0 if use_multiplier else 10.0 * multiplier
        row = [(0, 0, 1), 1.0, conductance]
        if kind == "RIV":
            row.append(1.0)
        row.append(0.0)
        aux = ["CONCENTRATION"]
        options = {}
        if use_multiplier:
            row.append(multiplier)
            aux.append("MULTIPLIER")
            options["auxmultname"] = "MULTIPLIER"
        constructor = {
            "GHB": flopy.mf6.ModflowGwfghb,
            "RIV": flopy.mf6.ModflowGwfriv,
            "DRN": flopy.mf6.ModflowGwfdrn,
        }[kind]
        constructor(gwf, pname=kind, auxiliary=aux, stress_period_data=[row], **options)
        # Default fully implicit BUY; this test does not depend on RHS iteration.
        flopy.mf6.ModflowGwfbuy(
            gwf,
            packagedata=[(0, 1.3, 0.0, transport, "CONCENTRATION")],
            density_filerecord=f"{flow}.density",
        )
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=f"{flow}.hds",
            budget_filerecord=f"{flow}.cbc",
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )
        gwt = flopy.mf6.ModflowGwt(sim, modelname=transport)
        imst = flopy.mf6.ModflowIms(
            sim,
            filename=f"{transport}.ims",
            outer_dvclose=1e-10,
            inner_dvclose=1e-10,
            rcloserecord=1e-9,
            linear_acceleration="BICGSTAB",
        )
        sim.register_ims_package(imst, [transport])
        flopy.mf6.ModflowGwtdis(gwt, **grid)
        flopy.mf6.ModflowGwtic(gwt, strt=concentration)
        flopy.mf6.ModflowGwtmst(gwt, porosity=0.35)
        flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
        flopy.mf6.ModflowGwtssm(
            gwt,
            sources=[("CHD", "AUX", "CONCENTRATION"), (kind, "AUX", "CONCENTRATION")],
        )
        flopy.mf6.ModflowGwtoc(
            gwt,
            concentration_filerecord=f"{transport}.ucn",
            saverecord=[("CONCENTRATION", "ALL")],
        )
        flopy.mf6.ModflowGwfgwt(
            sim, exgmnamea=flow, exgmnameb=transport, filename=f"{suffix}.gwfgwt"
        )
    return sim, None


def read_result(ws, suffix, kind):
    flow, transport = f"flow_{suffix}", f"gwt_{suffix}"
    heads = (
        flopy.utils.HeadFile(ws / f"{flow}.hds", precision="double").get_data().ravel()
    )
    density = (
        flopy.utils.HeadFile(ws / f"{flow}.density", text="DENSITY", precision="double")
        .get_data()
        .ravel()
    )
    conc = (
        flopy.utils.HeadFile(
            ws / f"{transport}.ucn", text="CONCENTRATION", precision="double"
        )
        .get_data()
        .ravel()
    )
    budget = flopy.utils.CellBudgetFile(ws / f"{flow}.cbc", precision="double")
    qboundary = budget.get_data(text=kind)[0]["q"].item()
    qchd = budget.get_data(text="CHD")[0]["q"].item()
    # CSR row1 is [diagonal, neighbor0]; positive rate is into cell1.
    face = budget.get_data(text="FLOW-JA-FACE")[0].ravel()[3]
    return dict(
        heads=heads,
        density=density,
        concentration=conc,
        boundary_flux=qboundary,
        chd_flux=qchd,
        face_flux=face,
    )


def physical_errors(result, kind, multiplier, concentration):
    h, c = result["heads"], result["concentration"]
    rho = 1.0 + 1.3 * concentration / 1000.0
    rb = rho if kind == "DRN" else 1.0
    cond = 10.0 * multiplier
    # Pressure at common z=.5. DRN effluent density equals aquifer density.
    expected = (2.0 * rho + cond * (rb + (rho - rb) * 0.5)) / (rho * (1.0 + cond))
    assert expected > 1.0, "River/drain must remain on its active branch"
    q = cond * (rb * (1.0 - 0.5) - rho * (h[1] - 0.5))
    face = rho * (h[0] - h[1])
    return dict(
        density_error=float(np.max(abs(result["density"] - 1000.0 * rho))),
        head_error=float(max(abs(h[0] - 2.0), abs(h[1] - expected))),
        pressure_flux_error=abs(q - result["boundary_flux"]),
        internal_flux_error=abs(face - result["face_flux"]),
        chd_flux_error=abs(face - result["chd_flux"]),
        cell_flow_residual=abs(face + q),
        water_budget_residual=abs(result["boundary_flux"] + result["chd_flux"]),
        solute_cell0_residual=abs(
            0.35 * (c[0] - concentration)
            - result["chd_flux"] * concentration
            + face * c[0]
        ),
        solute_cell1_residual=abs(
            0.35 * (c[1] - concentration) - face * c[0] - result["boundary_flux"] * c[1]
        ),
    )


def check_output(idx, test):
    kind, multiplier, concentration = settings[idx]
    results = [
        read_result(test.workspace, suffix, kind) for suffix in ["mult", "direct"]
    ]
    for result in results:
        errors = physical_errors(result, kind, multiplier, concentration)
        assert max(errors.values()) < 1e-8, errors
        if multiplier == 0.0:
            np.testing.assert_allclose(
                result["boundary_flux"], 0.0, rtol=0.0, atol=1e-10
            )
    for key in results[0]:
        np.testing.assert_allclose(
            results[0][key], results[1][key], rtol=0.0, atol=1e-8
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

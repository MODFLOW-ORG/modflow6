"""
Tests that the Buoyancy (BUY) and Viscosity (VSC) Packages are rejected when the
Variably Saturated Flow (VSF) Package is active in the same GWF model.  VSF
(Richards flow) is not supported with density- or viscosity-dependent flow, so
the run must terminate with a clear error rather than producing undefined
results.
"""

import os

import flopy
import pytest
from framework import TestFramework
from gwf_test_utils import get_vsf_soil_data

cases = ["vsf-buy", "vsf-vsc"]
package = ["buy", "vsc"]


def build_models(idx, test):
    nlay, nrow, ncol = 3, 1, 1
    delr = delc = 10.0
    column_height = 30.0
    delz = column_height / nlay
    top = delz
    botm = [top - (ilay + 1) * delz for ilay in range(nlay)]
    hk = 0.00944

    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim)

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="gwf", save_flows=True, newtonoptions="newton"
    )
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB"), ["gwf"]
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5 * (top + botm[-1]))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-4, sy=0.3, iconvert=0, steady_state=True)

    soil_data = get_vsf_soil_data("Celia1990-eq10-Haverkamp")
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        soil_model="Haverkamp",
        porosity=soil_data["porosity"],
        satres=soil_data["satres"],
        alphahvk=soil_data["alpha"],
        nhvk=soil_data["n"],
        betahvk=soil_data["beta"],
        khvk=soil_data["k"],
    )
    # carry an auxiliary concentration so BUY/VSC can reference it
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary="CONCENTRATION",
        stress_period_data=[[(0, 0, 0), top, 0.0], [(nlay - 1, 0, 0), botm[-1], 0.0]],
    )

    if package[idx] == "buy":
        flopy.mf6.ModflowGwfbuy(
            gwf, packagedata=[(0, 0.7, 0.0, "none", "CONCENTRATION")]
        )
    else:
        flopy.mf6.ModflowGwfvsc(
            gwf,
            viscref=1.0,
            nviscspecies=1,
            packagedata=[(0, 0.0, 20.0, "none", "CONCENTRATION")],
        )
    return sim


def check_output(idx, test):
    # the run must terminate because BUY/VSC are not supported with VSF
    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        normalized = " ".join(f.read().split())
    name = "Buoyancy (BUY)" if package[idx] == "buy" else "Viscosity (VSC)"
    tag = f"The {name} Package cannot be used with the Variably"
    assert tag in normalized, f"expected {name}-not-supported-with-VSF error not found"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=True,
    )
    test.run()

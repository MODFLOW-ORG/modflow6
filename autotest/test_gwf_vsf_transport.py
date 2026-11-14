"""
Tests that the Variably Saturated Flow (VSF) Package is rejected when the GWF
model it belongs to is coupled to a transport model.  VSF (Richards flow) is not
supported in combination with GWT or GWE, so the GWF-GWT and GWF-GWE exchanges
must terminate the run with a clear error rather than producing undefined
results.
"""

import os

import flopy
import pytest
from framework import TestFramework
from gwf_test_utils import get_vsf_soil_data

cases = ["vsf-gwt", "vsf-gwe"]
transport = ["gwt", "gwe"]


def build_gwf(sim, gwfname, nlay, nrow, ncol, delr, delc, top, botm, hk):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, newtonoptions="newton"
    )
    imsgwf = flopy.mf6.ModflowIms(
        sim, linear_acceleration="BICGSTAB", filename=f"{gwfname}.ims"
    )
    sim.register_ims_package(imsgwf, [gwfname])

    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5 * (top + botm[-1]))
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-4, sy=0.3, iconvert=0, steady_state=True)

    # variably saturated (Richards) flow
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
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), top], [(nlay - 1, 0, 0), botm[-1]]]
    )
    return gwf


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

    gwfname = "gwf"
    build_gwf(sim, gwfname, nlay, nrow, ncol, delr, delc, top, botm, hk)

    if transport[idx] == "gwt":
        gwtname = "gwt"
        gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)
        sim.register_ims_package(
            flopy.mf6.ModflowIms(sim, filename=f"{gwtname}.ims"), [gwtname]
        )
        flopy.mf6.ModflowGwtdis(
            gwt,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGwtic(gwt, strt=0.0)
        flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
        flopy.mf6.ModflowGwtssm(gwt)
        flopy.mf6.ModflowGwfgwt(
            sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
        )
    else:
        gwename = "gwe"
        gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename)
        sim.register_ims_package(
            flopy.mf6.ModflowIms(sim, filename=f"{gwename}.ims"), [gwename]
        )
        flopy.mf6.ModflowGwedis(
            gwe,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        flopy.mf6.ModflowGweic(gwe, strt=0.0)
        flopy.mf6.ModflowGweest(
            gwe,
            porosity=0.3,
            heat_capacity_water=4183.0,
            density_water=1000.0,
            latent_heat_vaporization=2453500.0,
            heat_capacity_solid=700.0,
            density_solid=2650.0,
        )
        flopy.mf6.ModflowGwessm(gwe)
        flopy.mf6.ModflowGwfgwe(
            sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
        )
    return sim


def check_output(idx, test):
    # the run must terminate because VSF is not supported with transport
    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        normalized = " ".join(f.read().split())
    tag = "Variably Saturated Flow (VSF) Package, which is not supported"
    assert tag in normalized, (
        "expected VSF-not-supported-with-transport error not found in mfsim.lst"
    )


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

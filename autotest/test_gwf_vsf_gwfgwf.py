"""
Tests that a GWF-GWF exchange is rejected when either connected GWF model uses
the Variably Saturated Flow (VSF) Package.  VSF (Richards flow) is not yet
supported across a GWF-GWF exchange (including LGR), so the run must terminate
with a clear error rather than producing undefined results.
"""

import os

import flopy
import pytest
from framework import TestFramework
from gwf_test_utils import get_vsf_soil_data

cases = ["vsf-gwfgwf"]


def add_vsf_model(sim, modelname, top, botm):
    hk = 0.00944
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=modelname, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5 * (top + botm))
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
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), top]])
    return gwf


def add_plain_model(sim, modelname, top, botm):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=modelname, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=botm)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), botm]])
    return gwf


def build_models(idx, test):
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim)
    # both models share a single solution, as required by a GWF-GWF exchange
    sim.register_ims_package(
        flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB"),
        ["gwf1", "gwf2"],
    )

    # gwf1 (VSF) stacked directly above gwf2 (plain)
    add_vsf_model(sim, "gwf1", top=0.0, botm=-1.0)
    add_plain_model(sim, "gwf2", top=-1.0, botm=-2.0)

    # vertical connection between the two single-cell models
    exchangedata = [[(0, 0, 0), (0, 0, 0), 0, 0.5, 0.5, 1.0, 0.0, 1.0]]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exchangedata),
        exgmnamea="gwf1",
        exgmnameb="gwf2",
        exchangedata=exchangedata,
        auxiliary=["ANGLDEGX", "CDIST"],
    )
    return sim


def check_output(idx, test):
    # the run must terminate because VSF is not supported across an exchange
    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        normalized = " ".join(f.read().split())
    tag = "Variably Saturated Flow (VSF) Package is inactive in both"
    assert tag in normalized, (
        "expected VSF-not-supported-with-exchange error not found in mfsim.lst"
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

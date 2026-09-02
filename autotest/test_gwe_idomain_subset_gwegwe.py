"""
Tests a GWE-GWE exchange in which one transport model simulates transport on a
subset of the cells of the flow model that supplies its flows.

Two flow models, split at the middle of a two-row grid, are coupled by a GWF-GWF
exchange.  Two pairs of energy transport models are coupled to them: a
full-domain pair that provides the reference solution, and a pair whose
downgradient model excludes the last thirty columns.  Excluding those cells
shifts the reduced node numbers of the exchange cells, so the flow and transport
exchanges match only when their cells are compared as user node numbers.  The
CND package is given no conductivity so that the energy transport is purely
advective and the reduced domain reproduces the full-domain solution exactly.

Cases:
  - gwegweexg : the flow models are coupled by a classic GWF-GWF exchange.
  - gwegweifm : the flow models are coupled through an interface model.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwegweexg", "gwegweifm"]
ifmodel = [False, True]

nlay, nrow, ncol = 1, 2, 50
delr, delc = 1.0, 1.0
top, botm = 1.0, [0.0]
strt = 1.0
hk = 1.0
porosity = 0.1
cpw, rhow = 4184.0, 1000.0
cps, rhos = 800.0, 2650.0

# the thermal front is retarded by about a factor of six relative to the water,
# so the simulation is run long enough for it to reach the transport boundary
perlen, nstp = 40.0, 400

# last active column of the downgradient transport model
jsub = 19

# one exchange connection per row, between the last column of the upgradient
# model and the first column of the downgradient model
exgdata = [
    [(0, i, ncol - 1), (0, i, 0), 1, 0.5 * delr, 0.5 * delr, delc, 0.0, delr]
    for i in range(nrow)
]


def sub_idomain():
    idomain = np.zeros((nlay, nrow, ncol), dtype=int)
    idomain[0, :, : jsub + 1] = 1
    return idomain


def add_dis(model, dis, xorigin):
    dis(
        model,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xorigin,
    )


def build_flow_model(sim, name, xorigin, upgradient):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, model_nam_file=f"{name}.nam"
    )
    add_dis(gwf, flopy.mf6.ModflowGwfdis, xorigin)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=hk, k33=hk, save_specific_discharge=True
    )
    if upgradient:
        flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data={0: [[(0, i, 0), 1.0, 1.0] for i in range(nrow)]},
            auxiliary="TEMPERATURE",
            pname="WEL-1",
        )
    else:
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: [[(0, i, ncol - 1), 0.0] for i in range(nrow)]},
            pname="CHD-1",
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def build_transport_model(sim, name, xorigin, idomain, upgradient):
    gwe = flopy.mf6.MFModel(
        sim, model_type="gwe6", modelname=name, model_nam_file=f"{name}.nam"
    )
    gwe.name_file.save_flows = True
    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xorigin,
        idomain=idomain,
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="upstream")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=porosity,
        heat_capacity_water=cpw,
        density_water=rhow,
        heat_capacity_solid=cps,
        density_solid=rhos,
    )
    # a CND package with no conductivity leaves the energy transport purely
    # advective; it is present because the GWE-GWE interface model takes its
    # equation scaling factor from CND
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.0, kts=0.0)
    flopy.mf6.ModflowGwefmi(gwe, flow_imbalance_correction=True)
    if upgradient:
        sources = [("WEL-1", "AUX", "TEMPERATURE")]
    else:
        sources = None
    flopy.mf6.ModflowGwessm(gwe, sources=sources)
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{name}.cbc",
        temperature_filerecord=f"{name}.ucn",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "LAST")],
    )
    return gwe


def add_transport_pair(sim, prefix, idomain2, gwfnames, hclose, rclose, nouter, ninner):
    names = [prefix + "a", prefix + "b"]
    build_transport_model(sim, names[0], 0.0, 1, True)
    build_transport_model(sim, names[1], ncol * delr, idomain2, False)
    flopy.mf6.ModflowGwegwe(
        sim,
        exgtype="GWE6-GWE6",
        gwfmodelname1=gwfnames[0],
        gwfmodelname2=gwfnames[1],
        nexg=len(exgdata),
        exgmnamea=names[0],
        exgmnameb=names[1],
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename=f"{prefix}.gwegwe",
    )
    for gwfname, gwename in zip(gwfnames, names):
        flopy.mf6.ModflowGwfgwe(
            sim,
            exgtype="GWF6-GWE6",
            exgmnamea=gwfname,
            exgmnameb=gwename,
            filename=f"{gwename}.gwfgwe",
        )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename=f"{prefix}.ims",
    )
    sim.register_ims_package(ims, names)
    return names


def build_models(idx, test):
    name = cases[idx]
    nouter, ninner = 100, 300
    hclose, rclose = 1e-8, 1e-8

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )

    gwfnames = ["gwfa", "gwfb"]
    build_flow_model(sim, gwfnames[0], 0.0, True)
    build_flow_model(sim, gwfnames[1], ncol * delr, False)
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata),
        exgmnamea=gwfnames[0],
        exgmnameb=gwfnames[1],
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwfgwf.gwfgwf",
        dev_interfacemodel_on=ifmodel[idx],
    )
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename="gwf.ims",
    )
    sim.register_ims_package(imsgwf, gwfnames)

    # full-domain transport pair used as the reference solution
    add_transport_pair(sim, "full", 1, gwfnames, hclose, rclose, nouter, ninner)
    add_transport_pair(
        sim, "sub", sub_idomain(), gwfnames, hclose, rclose, nouter, ninner
    )

    return sim, None


def temperature(test, name):
    fpth = os.path.join(test.workspace, f"{name}.ucn")
    tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
    return tobj.get_data()


def check_output(idx, test):
    ca_sub, cb_sub = temperature(test, "suba"), temperature(test, "subb")
    ca_full, cb_full = temperature(test, "fulla"), temperature(test, "fullb")

    # the upgradient model is unaffected by the reduced downgradient domain
    assert np.allclose(ca_sub, ca_full, atol=1e-9), (
        "upgradient temperatures do not match the full-domain solution"
    )

    # inactive transport cells are written as hnoflo
    inactive = sub_idomain() == 0
    assert np.all(cb_sub[inactive] == 1e30), "inactive cells were not written as hnoflo"

    # the retained cells of the downgradient model must match as well
    assert np.allclose(
        cb_sub[0, :, : jsub + 1], cb_full[0, :, : jsub + 1], atol=1e-9
    ), "downgradient temperatures do not match the full-domain solution"

    # the plume must have reached the truncation boundary
    assert cb_sub[0, 0, jsub] > 0.1, "plume did not reach the transport boundary"


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

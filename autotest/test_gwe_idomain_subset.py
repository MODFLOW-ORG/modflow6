"""
Tests a GWE model that simulates transport on a subset of the GWF model cells.

A one-dimensional flow model with an injection well at the upgradient end and a
constant head at the downgradient end supplies flows to one or two energy
transport models through separate GWF-GWE exchanges.  Conduction is omitted so
that the reduced domain reproduces the full-domain solution exactly.

Cases:
  - idomsub : the transport domain drops the downgradient cells; temperatures
              must match a full-domain transport model in the retained cells and
              the transport budget must close.
  - idomint : the transport domain drops cells at both ends and holds the
              upgradient cell at a specified temperature; temperatures must
              again match the full-domain model in the retained cells.
  - idomnot : a cell that is active in the transport model and inactive in the
              flow model is an error (xfail).
  - idomapt : a transport domain that excludes a cell connected to an advanced
              package handled by MWE is an error (xfail).
  - idommwe : the transport domain drops the upgradient cells and keeps a
              well handled by MWE with well-bore conduction, so the transport
              and flow node numbers of the well cell differ; temperatures must
              match the full-domain model in the retained cells.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gweidomsub", "gweidomint", "gweidomnot", "gweidomapt", "gweidommwe"]

nlay, nrow, ncol = 1, 1, 100
delr, delc = 1.0, 1.0
top, botm = 1.0, [0.0]
strt = 1.0
hk = 1.0
ss = 1.0e-5
porosity = 0.1
cpw, rhow = 4184.0, 1000.0
cps, rhos = 800.0, 2650.0

# the thermal front is retarded by about a factor of six relative to the water,
# so the simulation is run long enough for it to move past the transport domain
perlen, nstp = 50.0, 500

# first and last active transport column for each case
icol0 = [0, 20, 0, 0, 20]
icol1 = [59, 59, 59, 59, 59]

# column deactivated in the flow model, or None
icolgwf = [None, None, 80, None, None]

# column with a multi-aquifer well outside the transport domain, or None
icolmaw = [None, None, None, 80, 40]

# multi-aquifer well injection rate
mawrate = [0.0, 0.0, 0.0, 0.0, 0.5]

# column with a specified transport temperature, or None
icolctp = [None, 20, None, None, 20]


def transport_idomain(idx):
    idomain = np.zeros((nlay, nrow, ncol), dtype=int)
    idomain[0, 0, icol0[idx] : icol1[idx] + 1] = 1
    if icolgwf[idx] is not None:
        idomain[0, 0, icolgwf[idx]] = 1
    return idomain


def build_flow_model(sim, gwfname, idx, hclose, rclose, nouter, ninner):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, save_flows=True, model_nam_file=f"{gwfname}.nam"
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwf.name])

    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    if icolgwf[idx] is not None:
        idomain[0, 0, icolgwf[idx]] = 0
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=hk, k33=hk, save_specific_discharge=True
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=0, ss=ss, sy=0.0, steady_state=None, transient={0: True}
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, ncol - 1), 0.0]]},
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={0: [[(0, 0, 0), 1.0, 1.0]]},
        auxiliary="TEMPERATURE",
        pname="WEL-1",
    )
    if icolmaw[idx] is not None:
        flopy.mf6.ModflowGwfmaw(
            gwf,
            nmawwells=1,
            packagedata=[[0, 0.1, 0.0, 1.0, "THIEM", 1]],
            connectiondata=[[0, 0, (0, 0, icolmaw[idx]), 1.0, 0.0, 1.0, 0.1]],
            perioddata={0: [[0, "RATE", mawrate[idx]]]},
            pname="MAW-1",
        )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def build_transport_model(
    sim, gwename, idomain, idx, hclose, rclose, nouter, ninner, imaw=False
):
    gwe = flopy.mf6.MFModel(
        sim,
        model_type="gwe6",
        modelname=gwename,
        model_nam_file=f"{gwename}.nam",
    )
    gwe.name_file.save_flows = True
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(ims, [gwe.name])

    flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
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
    flopy.mf6.ModflowGwefmi(gwe, flow_imbalance_correction=True)
    flopy.mf6.ModflowGwessm(gwe, sources=[("WEL-1", "AUX", "TEMPERATURE")])
    if icolctp[idx] is not None:
        flopy.mf6.ModflowGwectp(
            gwe,
            stress_period_data={0: [[(0, 0, icolctp[idx]), 1.0]]},
            pname="CTP-1",
        )
    if imaw:
        flopy.mf6.ModflowGwemwe(
            gwe,
            flow_package_name="MAW-1",
            packagedata=[[0, 0.0, 1.0, 1.0]],
            mweperioddata={0: [[0, "RATE", 1.0]]},
            pname="MWE-1",
        )
    flopy.mf6.ModflowGweoc(
        gwe,
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )
    return gwe


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

    gwfname = "gwf_" + name
    build_flow_model(sim, gwfname, idx, hclose, rclose, nouter, ninner)

    # full-domain transport model used as the reference solution
    if idx in (0, 1, 4):
        gwename = "full_" + name
        build_transport_model(
            sim,
            gwename,
            np.ones((nlay, nrow, ncol), dtype=int),
            idx,
            hclose,
            rclose,
            nouter,
            ninner,
            imaw=icolmaw[idx] is not None,
        )
        flopy.mf6.ModflowGwfgwe(
            sim,
            exgtype="GWF6-GWE6",
            exgmnamea=gwfname,
            exgmnameb=gwename,
            filename=f"{gwename}.gwfgwe",
        )

    gwename = "sub_" + name
    build_transport_model(
        sim,
        gwename,
        transport_idomain(idx),
        idx,
        hclose,
        rclose,
        nouter,
        ninner,
        imaw=icolmaw[idx] is not None,
    )
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{gwename}.gwfgwe",
    )

    return sim, None


def budget_discrepancy(test, gwename):
    """Largest relative energy budget discrepancy in the model listing file."""
    with open(os.path.join(test.workspace, f"{gwename}.lst")) as f:
        text = f.read()
    residual = np.array(re.findall(r"IN - OUT =\s*(\S+)", text), dtype=float)
    total = np.array(re.findall(r"TOTAL IN =\s*(\S+)", text), dtype=float)
    assert residual.size > 0, f"no budget written to {gwename}.lst"
    return np.abs(residual) / np.where(total == 0.0, 1.0, total)


def check_output(idx, test):
    name = cases[idx]

    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        text = " ".join(f.read().split())

    if idx == 2:
        srch = "is active in the GWE Model but is not active in the GWF Model"
        assert srch in text, (
            "expected an error for a GWE idomain that is not a subset of the "
            "GWF idomain"
        )
        return

    if idx == 3:
        srch = "which is not active in the GWE Model"
        assert srch in text, (
            "expected an error for an advanced package cell that is outside "
            "the transport domain"
        )
        return

    gwename = "sub_" + name
    fpth = os.path.join(test.workspace, f"{gwename}.ucn")
    tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
    tsub = tobj.get_data()

    # the transport budget must close
    pd = budget_discrepancy(test, gwename)
    assert pd.max() < 1e-6, (
        f"transport budget does not close, max discrepancy {pd.max()}"
    )

    # inactive transport cells are written as hnoflo
    inactive = transport_idomain(idx) == 0
    assert np.all(tsub[inactive] == 1e30), "inactive cells were not written as hnoflo"

    # temperatures in the retained cells must match the full-domain model
    fpth = os.path.join(test.workspace, f"full_{name}.ucn")
    tobj = flopy.utils.HeadFile(fpth, precision="double", text="TEMPERATURE")
    tfull = tobj.get_data()

    j0, j1 = icol0[idx], icol1[idx] + 1
    assert np.allclose(tsub[0, 0, j0:j1], tfull[0, 0, j0:j1], atol=1e-9), (
        "temperatures in the reduced transport domain do not match the "
        "full-domain solution"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        xfail=idx in (2, 3),
    )
    test.run()

"""
Test the specific discharge calculation for an LGR-like simulation that has
a parent model and a child model.  The child model is inset into the parent
model, but they both have the same resolution, so it is essentially a simple
3D grid.  The child qx velocity should be the same as the qx velocity in
the parent grid.  The heads are also compared.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.lgrutil import Lgr
from framework import TestFramework

cases = ["npf04"]
namea = "a"
nameb = "b"


def build_models(idx, test):
    # grid properties
    nlay = 3
    nrow = 6
    ncol = 6
    delr = 100.0
    delc = 100.0
    top = 300.0
    botm = [200.0, 100.0, 0.0]

    # hydraulic properties
    hk = 1.0
    vk = 1.0

    # Set the idomain of the parent model in order to
    # define where the child model will be located
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[:, 2:4, 2:4] = 0

    ncpp = 1
    ncppl = [1, 1, 1]
    lgr = Lgr(nlay, nrow, ncol, delr, delc, top, botm, idomain, ncpp, ncppl)

    name = cases[idx]

    # build MODFLOW 6 files
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)])

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=namea, save_flows=True)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(sim, outer_dvclose=1e-9, inner_dvclose=1.0e-9)

    # dis
    dis = flopy.mf6.ModflowGwfdis(
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

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=top)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_specific_discharge=True, icelltype=0, k=hk, k33=vk
    )

    # chd
    chdspd = []
    for k in range(nlay):
        for i in range(nrow):
            chdspd.append([(k, i, 0), 1.0])
            chdspd.append([(k, i, ncol - 1), 6.0])
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        pname="oc",
        budget_filerecord=f"{namea}.cbc",
        head_filerecord=f"{namea}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create child gwf model
    cgwf = flopy.mf6.ModflowGwf(sim, modelname=nameb, save_flows=True)
    cnlay, cnrow, cncol = lgr.get_shape()
    cdelr, cdelc = lgr.get_delr_delc()
    ctop, cbotm = lgr.get_top_botm()
    xorigin, yorigin = lgr.get_lower_left()
    cidomain = lgr.get_idomain()
    cdis = flopy.mf6.ModflowGwfdis(
        cgwf,
        nlay=cnlay,
        nrow=cnrow,
        ncol=cncol,
        delr=cdelr,
        delc=cdelc,
        top=ctop,
        botm=cbotm,
        idomain=cidomain,
        xorigin=xorigin,
        yorigin=yorigin,
    )
    cic = flopy.mf6.ModflowGwfic(cgwf, pname="ic", strt=top)
    cnpf = flopy.mf6.ModflowGwfnpf(
        cgwf, save_specific_discharge=True, icelltype=0, k=hk, k33=vk
    )
    oc = flopy.mf6.ModflowGwfoc(
        cgwf,
        pname="oc",
        budget_filerecord=f"{nameb}.cbc",
        head_filerecord=f"{nameb}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # exchange information
    exchangedata = lgr.get_exchange_data(angldegx=True, cdist=True)
    nexg = len(exchangedata)
    gwfe = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="gwf6-gwf6",
        exgmnamea="a",
        exgmnameb="b",
        nexg=nexg,
        auxiliary=[("angldegx", "cdist")],
        exchangedata=exchangedata,
    )

    return sim, None


def qxqyqz(fname, nlay, nrow, ncol):
    nodes = nlay * nrow * ncol
    cbb = flopy.utils.CellBudgetFile(fname, precision="double")
    spdis = cbb.get_data(text="DATA-SPDIS")[0]
    qx = np.ones((nodes), dtype=float) * 1.0e30
    qy = np.ones((nodes), dtype=float) * 1.0e30
    qz = np.ones((nodes), dtype=float) * 1.0e30
    n0 = spdis["node"] - 1
    qx[n0] = spdis["qx"]
    qy[n0] = spdis["qy"]
    qz[n0] = spdis["qz"]
    qx = qx.reshape(nlay, nrow, ncol)
    qy = qy.reshape(nlay, nrow, ncol)
    qz = qz.reshape(nlay, nrow, ncol)
    qx = np.ma.masked_equal(qx, 1.0e30)
    qy = np.ma.masked_equal(qy, 1.0e30)
    qz = np.ma.masked_equal(qz, 1.0e30)
    return qx, qy, qz


def check_output(idx, test):
    # make sure parent head is same as child head in same column
    fname = os.path.join(test.workspace, f"{namea}.hds")
    hdobj = flopy.utils.HeadFile(fname)
    ha = hdobj.get_data()
    fname = os.path.join(test.workspace, f"{nameb}.hds")
    hdobj = flopy.utils.HeadFile(fname)
    hb = hdobj.get_data()
    msg = f"Heads should be the same {ha[0, 1, 2]} {hb[0, 0, 0]}"
    assert np.allclose(ha[0, 1, 2], hb[0, 0, 0]), msg

    # make sure specific discharge is calculated correctly for child and
    # parent models (even though child model has same resolution as parent
    fname = os.path.join(test.workspace, f"{namea}.cbc")
    nlaya, nrowa, ncola = ha.shape
    qxa, qya, qza = qxqyqz(fname, nlaya, nrowa, ncola)
    fname = os.path.join(test.workspace, f"{nameb}.cbc")
    nlayb, nrowb, ncolb = hb.shape
    qxb, qyb, qzb = qxqyqz(fname, nlayb, nrowb, ncolb)
    msg = f"qx should be the same {qxa[0, 2, 1]} {qxb[0, 0, 0]}"
    assert np.allclose(qxa[0, 2, 1], qxb[0, 0, 0]), msg

    cbcpth = os.path.join(test.workspace, f"{namea}.cbc")
    grdpth = os.path.join(test.workspace, f"{namea}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(grdpth)
    cbb = flopy.utils.CellBudgetFile(cbcpth, precision="double")
    flow_ja_face = cbb.get_data(text="FLOW-JA-FACE")
    ia = grb._datadict["IA"] - 1
    for fjf in flow_ja_face:
        fjf = fjf.flatten()
        res = fjf[ia[:-1]]
        errmsg = f"min or max residual too large {res.min()} {res.max()}"
        assert np.allclose(res, 0.0, atol=1.0e-6), errmsg


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

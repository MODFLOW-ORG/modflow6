"""
Test the interface model approach, when running
with a GWF-GWF exchange and XT3D applied on it.
It compares the result for a simple LGR configuration
to the analytical values:

        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 0 0 0 1 1
(H=1.0) 1 1 0 0 0 1 1 (H=0.0)
        1 1 0 0 0 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1

with the region with ibound == 0 being simulated on the
a refined, 9x9 grid.

This is also the first test problem presented in
the MODFLOW-USG manual: 'test006_2models'

When running without XT3D, the results will disagree
with theory because the CVFD requirements are violated at the
at the LGR interface. We compare heads, specific discharge, and
confirm that there is no budget error.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.lgrutil import Lgr
from framework import TestFramework

cases = ["ifmod_xt3d01"]

# globally for convenience...
useXT3D = True
parent_name = "parent"
child_name = "child"
h_left = 1.0
h_right = 0.0
delr = 100.0
delc = 100.0
k11 = 1.0
k33 = 1.0
child_domain = None
hclose = None


def get_model(idx, dir):
    global child_domain
    global hclose

    name = cases[idx]

    # tdis period data
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-12, 1e-3, 0.97
    h_start = 1.0

    # dis
    nlay, nrow, ncol = 1, 7, 7

    row_s, row_e = 3, 5
    col_s, col_e = 3, 5

    ref_fct = 5
    nrowc = ref_fct * ((row_e - row_s) + 1)
    ncolc = ref_fct * ((col_e - col_s) + 1)

    idomain = np.ones((nlay, nrow, ncol))
    idomain[:, row_s - 1 : row_e, col_s - 1 : col_e] = 0

    delrc = delr / ref_fct
    delcc = delc / ref_fct
    tops = [0.0, -100.0]

    xoriginc = 2 * delr
    yoriginc = 2 * delc

    xmin = 0.0
    xmax = ncol * delr
    ymin = 0.0
    ymax = nrow * delc
    model_domain = [xmin, xmax, ymin, ymax]
    xminc = xoriginc
    xmaxc = xoriginc + ncolc * delrc
    yminc = yoriginc
    ymaxc = yoriginc + nrowc * delcc
    child_domain = [xminc, xmaxc, yminc, ymaxc]

    # boundary stress period data
    left_chd = [
        [(ilay, irow, 0), h_left] for ilay in range(nlay) for irow in range(nrow)
    ]
    right_chd = [
        [(ilay, irow, ncol - 1), h_right]
        for ilay in range(nlay)
        for irow in range(nrow)
    ]
    chd_data = left_chd + right_chd
    chd_spd = {0: chd_data}

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="ALL",
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )

    # The parent model:
    gwf = flopy.mf6.ModflowGwf(sim, modelname=parent_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
        idomain=idomain,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=True,
        save_flows=True,
        icelltype=0,
        k=k11,
        k33=k33,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{parent_name}.hds",
        budget_filerecord=f"{parent_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # The child model:
    gwfc = flopy.mf6.ModflowGwf(sim, modelname=child_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwfc,
        nlay=nlay,
        nrow=nrowc,
        ncol=ncolc,
        delr=delrc,
        delc=delcc,
        top=tops[0],
        botm=tops[1:],
        xorigin=xoriginc,
        yorigin=yoriginc,
    )
    ic = flopy.mf6.ModflowGwfic(gwfc, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwfc,
        save_specific_discharge=True,
        xt3doptions=True,
        save_flows=True,
        icelltype=0,
        k=k11,
        k33=k33,
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwfc,
        head_filerecord=f"{child_name}.hds",
        budget_filerecord=f"{child_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # LGR:
    nrowp = gwf.dis.nrow.get_data()
    ncolp = gwf.dis.ncol.get_data()
    delrp = gwf.dis.delr.array
    delcp = gwf.dis.delc.array
    topp = gwf.dis.top.array
    botmp = gwf.dis.botm.array
    idomainp = gwf.dis.idomain.array

    lgr = Lgr(
        nlay, nrowp, ncolp, delrp, delcp, topp, botmp, idomainp, ncpp=ref_fct, ncppl=1
    )

    exgdata = lgr.get_exchange_data(angldegx=True, cdist=True)
    exgdata_withbname = []
    for exg in exgdata:
        l = exg
        angle = l[-2]
        if angle == 0:
            bname = "left"
        elif angle == 90.0:
            bname = "bottom"
        elif angle == 180.0:
            bname = "right"
        elif angle == 270.0:
            bname = "top"
        l.append(bname)
        exgdata_withbname.append(l)

    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata),
        exgmnamea=parent_name,
        exgmnameb=child_name,
        exchangedata=exgdata,
        xt3d=useXT3D,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=True,
        auxiliary=["ANGLDEGX", "CDIST"],
    )
    obslist = []
    obstype = "FLOW-JA-FACE"
    for iexg, exg_connection in enumerate(exgdata):
        obslist.append((f"f{iexg + 1}", obstype, (iexg,)))

    gwfgwfobs = {}
    gwfgwfobs["gwf_obs.csv"] = obslist
    gwfgwfobs["gwf_obs_boundnames.csv"] = [
        ["OBSLEFT", "FLOW-JA-FACE", "LEFT"],
        ["OBSRIGHT", "FLOW-JA-FACE", "RIGHT"],
        ["OBSTOP", "FLOW-JA-FACE", "TOP"],
        ["OBSBOTTOM", "FLOW-JA-FACE", "BOTTOM"],
    ]

    fname = "gwfgwf.obs"
    gwfgwf.obs.initialize(
        filename=fname, digits=25, print_input=True, continuous=gwfgwfobs
    )

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
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
    print("comparing heads and spec. discharges to analytical result...")

    fpth = os.path.join(test.workspace, f"{parent_name}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()

    fpth = os.path.join(test.workspace, f"{parent_name}.cbc")
    nlay, nrow, ncol = heads.shape
    qxb, qyb, qzb = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{child_name}.hds")
    hds_c = flopy.utils.HeadFile(fpth)
    heads_c = hds_c.get_data()

    fpth = os.path.join(test.workspace, f"{child_name}.cbc")
    nlay, nrow, ncol = heads_c.shape
    qxb_c, qyb_c, qzb_c = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{parent_name}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    mg = grb.modelgrid

    fpth = os.path.join(test.workspace, f"{child_name}.dis.grb")
    grb_c = flopy.mf6.utils.MfGrdFile(fpth)
    mg_c = grb_c.modelgrid

    xyc = mg.xycenters
    xyc_c = mg_c.xycenters

    # the exact results:
    xleft = xyc[0][0]
    xright = xyc[0][-1]

    def exact(x):
        return h_left + (h_right - h_left) * (x - xleft) / (xright - xleft)

    qx_exact = -k11 * (h_right - h_left) / ((mg.ncol - 1) * delr)

    # first check the parent
    for irow in range(mg.nrow):
        for icol in range(mg.ncol):
            xc = xyc[0][icol]
            h = heads[0, irow, icol]
            if h != 1e30:
                diff = abs(h - exact(xc))
                assert diff < 10 * hclose, (
                    f"head difference in parent model {diff}"
                    f" exceeds solver tolerance (x10) {10 * hclose}"
                    f" for row {irow + 1} and col {icol + 1}\n"
                    f"(should be {exact(xc)}, was {h})"
                )

    for irow in range(mg.nrow):
        for icol in range(mg.ncol):
            if qyb[0, irow, icol] != "--":
                diff = abs(qyb[0, irow, icol])
                assert diff < 10 * hclose, (
                    "Specific discharge should not have a y-component in this model"
                )
            if qxb[0, irow, icol] != "--":
                diff = abs(qxb[0, irow, icol] - qx_exact)
                assert diff < 10 * hclose, (
                    f"Difference in spec. dis. for parent {diff}"
                    f" exceeds solver tolerance (x10) {10 * hclose}"
                    f" for row {irow + 1} and col {icol + 1}"
                )

    # and now the child
    for irow in range(mg_c.nrow):
        for icol in range(mg_c.ncol):
            xc = xyc_c[0][icol] + child_domain[0]
            h = heads_c[0, irow, icol]
            if h != 1e30:
                diff = abs(h - exact(xc))
                assert diff < 10 * hclose, (
                    f"Head difference in child model {diff}"
                    f" exceeds solver tolerance (x10) {10 * hclose}"
                    f" for row {irow + 1} and col {icol + 1}"
                )

    for irow in range(mg_c.nrow):
        for icol in range(mg_c.ncol):
            if qyb_c[0, irow, icol] != "--":
                diff = abs(qyb_c[0, irow, icol])
                assert diff < 10 * hclose, (
                    "Specific discharge should not have a y-component in this model"
                )
            if qxb_c[0, irow, icol] != "--":
                diff = abs(qxb_c[0, irow, icol] - qx_exact)
                assert diff < 10 * hclose, (
                    f"Difference in spec. dis. for child {diff}"
                    f" exceeds solver tolerance (x10) {10 * hclose}"
                    f" for row {irow + 1} and col {icol + 1}"
                )

    # todo: mflistbudget
    # check cumulative balance error from .lst file
    for mname in [parent_name, child_name]:
        fpth = os.path.join(test.workspace, f"{mname}.lst")
        for line in open(fpth):
            if line.lstrip().startswith("PERCENT"):
                cumul_balance_error = float(line.split()[3])
                assert abs(cumul_balance_error) < 0.00001, (
                    f"Cumulative balance error = {cumul_balance_error} for {mname}, "
                    "should equal 0.0"
                )

    # Check on residual, which is stored in diagonal position of
    # flow-ja-face.  Residual should be less than convergence tolerance,
    # or this means the residual term is not added correctly.
    fpth = os.path.join(test.workspace, f"{parent_name}.cbc")
    cbb = flopy.utils.CellBudgetFile(fpth)
    flow_ja_face = cbb.get_data(idx=0)
    assert len(flow_ja_face) > 0, (
        "Could not check residuals as flow-ja-face could not be found"
    )
    ia = grb._datadict["IA"] - 1
    for fjf in flow_ja_face:
        fjf = fjf.flatten()
        res = fjf[ia[:-1]]
        errmsg = f"min or max residual too large {res.min()} {res.max()}"
        assert np.allclose(res, 0.0, atol=1.0e-6), errmsg

    # Read gwf-gwf observations values
    fpth = os.path.join(test.workspace, "gwf_obs.csv")
    with open(fpth) as f:
        lines = f.readlines()
    obsnames = [name for name in lines[0].strip().split(",")[1:]]
    obsvalues = [float(v) for v in lines[1].strip().split(",")[1:]]

    # Extract the gwf-gwf flows stored in parent budget file
    fpth = os.path.join(test.workspace, f"{parent_name}.cbc")
    cbb = flopy.utils.CellBudgetFile(fpth, precision="double")
    parent_exchange_flows = cbb.get_data(
        kstpkper=(0, 0), text="FLOW-JA-FACE", paknam="GWF-GWF_1"
    )[0]
    parent_exchange_flows = parent_exchange_flows["q"]

    # Extract the gwf-gwf flows stored in child budget file
    fpth = os.path.join(test.workspace, f"{child_name}.cbc")
    cbb = flopy.utils.CellBudgetFile(fpth, precision="double")
    child_exchange_flows = cbb.get_data(
        kstpkper=(0, 0), text="FLOW-JA-FACE", paknam="GWF-GWF_1"
    )[0]
    child_exchange_flows = child_exchange_flows["q"]

    # Ensure observations are the same as parent exchange flows
    # and negative child exchange flows
    assert np.allclose(obsvalues, parent_exchange_flows), (
        "exchange observations do not match parent exchange flows"
    )
    assert np.allclose(obsvalues, -child_exchange_flows), (
        "exchange observations do not match child exchange flows"
    )

    # Read the lumped boundname observations values
    fpth = os.path.join(test.workspace, "gwf_obs_boundnames.csv")
    with open(fpth) as f:
        lines = f.readlines()
    obsnames = [name for name in lines[0].strip().split(",")[1:]]
    obsvalues = [float(v) for v in lines[1].strip().split(",")[1:]]
    assert np.allclose(obsvalues, [-50.0, 50.0, 0, 0.0]), (
        "boundname observations do not match expected results"
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

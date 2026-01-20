"""
Test the interface model approach for multiple (2) exchanges between
the same two models. One exchange has XT3D and the other one doesn't.

        'parent: 1x6x8'          'child: 1x16x16'

     1  1  1  1  1  1  1  1                      XT3D on exg
     1  1  0  0  0  0  1  1         1, ..., 1
     1  1  0  0  0  0  1  1         1, ..., 1
----------------------------------------------------------------
     1  1  0  0  0  0  1  1         1, ..., 1
     1  1  0  0  0  0  1  1         1, ..., 1
     1  1  1  1  1  1  1  1                     no XT3D on exg


with the top half of the exchange connections being part of
exchange_north and the others of exchange_south. The former
will have the XT3D calculation enabled.

TODO: (how) will this affect accuracy?
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.lgrutil import Lgr
from framework import TestFramework

cases = ["ifmod_mult_exg"]
name_parent = "parent"
name_child = "child"
g_delr = 10.0
g_k11 = 10.0
g_hleft = 120.0
g_hright = 190.0
g_hclose = 10e-12


def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = g_hclose, 1e-3, 0.97

    # refinement factor
    rft = 4

    # model spatial discretization for parent
    nlay = 1
    ncol = 8
    nrow = 6

    # cell spacing
    delr = g_delr
    delc = 10.0

    area = delr * delc

    # top/bot of the aquifer
    tops = [0.0, -5.0]

    # hydraulic conductivity
    k11 = g_k11

    # boundary stress period data
    h_left = g_hleft
    h_right = g_hright

    # initial head
    h_start = 135.0

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=dir,
        memory_print_option="all",
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

    # boundary data
    left_chd = [
        [(ilay, irow, 0), h_left] for irow in range(nrow) for ilay in range(nlay)
    ]
    right_chd = [
        [(ilay, irow, ncol - 1), h_right]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    chd_data = left_chd + right_chd
    chd_spd = {0: chd_data}

    # idomain
    idomain = np.ones((nlay, nrow, ncol))
    idomain[:, 1:5, 2:6] = 0

    # --------------------------------------
    # parent model
    # --------------------------------------
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name_parent, save_flows=True)
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
        save_flows=True,
        icelltype=0,
        k=k11,
        xt3doptions=False,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name_parent}.hds",
        budget_filerecord=f"{name_parent}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # --------------------------------------
    # child model (gwfc)
    # --------------------------------------
    gwfc = flopy.mf6.ModflowGwf(sim, modelname=name_child, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwfc,
        nlay=nlay,
        nrow=rft * 4,
        ncol=rft * 4,
        delr=delr / rft,
        delc=delc / rft,
        top=tops[0],
        botm=tops[1:],
        xorigin=2 * delr,
        yorigin=delc,
        angrot=0.0,
    )
    ic = flopy.mf6.ModflowGwfic(gwfc, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwfc,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
        xt3doptions=False,
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwfc,
        head_filerecord=f"{name_child}.hds",
        budget_filerecord=f"{name_child}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # --------------------------------------
    # LGR
    # --------------------------------------
    nrowp = gwf.dis.nrow.get_data()
    ncolp = gwf.dis.ncol.get_data()
    delrp = gwf.dis.delr.array
    delcp = gwf.dis.delc.array
    topp = gwf.dis.top.array
    botmp = gwf.dis.botm.array
    idomainp = gwf.dis.idomain.array

    lgr = Lgr(
        nlay, nrowp, ncolp, delrp, delcp, topp, botmp, idomainp, ncpp=rft, ncppl=1
    )

    exgdata = lgr.get_exchange_data(angldegx=True, cdist=True)
    exgdata_north = [e for e in exgdata if (e[0])[1] < 3]  # northern three rows
    exgdata_south = [e for e in exgdata if (e[0])[1] > 2]  # southern three rows

    # north, has XT3D
    flopy.mf6.ModflowGwfgwf(
        sim,
        dev_interfacemodel_on=True,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata_north),
        exgmnamea=name_parent,
        exgmnameb=name_child,
        exchangedata=exgdata_north,
        xt3d=True,
        print_flows=True,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="north_xt3d.gwfgwf",
    )

    # south, no XT3D
    flopy.mf6.ModflowGwfgwf(
        sim,
        dev_interfacemodel_on=True,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata_south),
        exgmnamea=name_parent,
        exgmnameb=name_child,
        exchangedata=exgdata_south,
        xt3d=False,
        print_flows=True,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="south_noxt3d.gwfgwf",
    )

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, f"{name_parent}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()

    fpth = os.path.join(test.workspace, f"{name_child}.hds")
    hds_c = flopy.utils.HeadFile(fpth)
    heads_c = hds_c.get_data()

    fpth = os.path.join(test.workspace, f"{name_parent}.dis.grb")
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    mg = grb.modelgrid

    fpth = os.path.join(test.workspace, f"{name_child}.dis.grb")
    grb_c = flopy.mf6.utils.MfGrdFile(fpth)
    mg_c = grb_c.modelgrid

    xyc = mg.xycenters
    xyc_c = mg_c.xycenters

    # the exact results:
    xleft = xyc[0][0]
    xright = xyc[0][-1]

    def exact(x):
        return g_hleft + (g_hright - g_hleft) * (x - xleft) / (xright - xleft)

    # first compare the parent's second row to exact result
    maxdiff_parent_north = 0.0
    for icol in range(mg.ncol):
        xc = xyc[0][icol]
        h = heads[0, 1, icol]  # second row
        if h != 1e30:
            diff = abs(h - exact(xc))
            if diff > maxdiff_parent_north:
                maxdiff_parent_north = diff

    maxdiff_parent_south = 0.0
    for icol in range(mg.ncol):
        xc = xyc[0][icol]
        h = heads[0, 4, icol]  # fifth row
        if h != 1e30:
            diff = abs(h - exact(xc))
            if diff > maxdiff_parent_south:
                maxdiff_parent_south = diff

    assert maxdiff_parent_south > maxdiff_parent_north

    # maxdiff_child_north = 0.0
    # for icol in range(mg_c.ncol):
    #     xc = xyc_c[0][icol] + 2*g_delr # shift by xorigin
    #     h = heads_c[0, 1, icol]  # second row
    #     if h != 1e30:
    #         diff = abs(h - exact(xc))
    #         if diff > maxdiff_child_north:
    #             maxdiff_child_north = diff
    #
    # maxdiff_child_south = 0.0
    # for icol in range(mg_c.ncol):
    #     xc = xyc_c[0][icol] + 2*g_delr # shift by xorigin
    #     h = heads_c[0, 14, icol]  # fifteenth row
    #     if h != 1e30:
    #         diff = abs(h - exact(xc))
    #         if diff > maxdiff_child_south:
    #             maxdiff_child_south = diff

    # assert maxdiff_child_south > maxdiff_child_north


@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

"""
Test the interface model approach.
It compares the result of a single, strongly anisotropic model
with XT3D enabled to the equivalent case where the domain is
decomposed into 4 models connected with GWF-GWF exchanges all
having XT3D enabled. Note the location of the well W, in the
bottom right corner of model "tl" (and also in "ref" of course)

          'ref'                    'tl'            'tr'

                                1 1 1 1 1        1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1        1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1    +   1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1        1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 W        1 1 1 1 1
   1 1 1 1 W 1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1    VS        +                +
   1 1 1 1 1 1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1        1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1        1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1    +   1 1 1 1 1
                                1 1 1 1 1        1 1 1 1 1
                                1 1 1 1 1        1 1 1 1 1

                                  'bl'             'br'

The head values should always be identical. All models are
part of the same solution for convenience.
In addition, a check on the x,y,z components of specific discharge
is present. Finally, the budget error is checked.
"""

import os
from types import SimpleNamespace

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ifmod_xt3d03"]

hclose_check = 1e-9
max_inner_it = 300
useXT3D = True

# model spatial discretization
nlay = 1
ncol = 10
ncol_split = 5
nrow = 10
nrow_split = 5

# cell spacing
delr = 10.0
delc = 10.0
area = delr * delc

# shift (hor. and vert.)
shift_some_x = -20 * delr  # avoids overlap
shift_x = 5 * delr
shift_y = 5 * delc

# top/bot of the aquifer
tops = [0.0, -5.0]

# hydraulic conductivity
k11 = 10.0
k22 = 0.1
k_angle = 45.0

# boundary stress period data
h_left = -2.0
h_right = -2.0

# initial head
h_start = -2.0

# well
well_id = (0, 4, 4)
well_rate = -1.0


def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, max_inner_it
    hclose, rclose, relax = hclose_check, 1e-3, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
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

    # reference model
    dis_params = SimpleNamespace(
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        xorigin=shift_some_x,
        yorigin=0.0,
        tops=tops,
    )
    create_gwf_model(sim, "ref", dis_params)

    # top-left model
    dis_params.nrow = nrow_split
    dis_params.ncol = ncol_split
    dis_params.xorigin = 0.0
    dis_params.yorigin = shift_y
    create_gwf_model(sim, "tl", dis_params)

    # bottom-left model
    dis_params.xorigin = 0.0
    dis_params.yorigin = 0.0
    create_gwf_model(sim, "bl", dis_params)

    # top-right model
    dis_params.xorigin = shift_x
    dis_params.yorigin = shift_y
    create_gwf_model(sim, "tr", dis_params)

    # bottom-right model
    dis_params.xorigin = shift_x
    dis_params.yorigin = 0.0
    create_gwf_model(sim, "br", dis_params)

    # two types of exchange data:
    # tl-tr, bl-br (0 deg)
    lr_data = [
        [
            (0, irow, ncol_split - 1),
            (0, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            0.0,
            delr,
        ]
        for irow in range(nrow_split)
    ]
    # tl-bl, tr-br (270 deg)
    tb_data = [
        [
            (0, nrow_split - 1, icol),
            (0, 0, icol),
            1,
            delc / 2.0,
            delc / 2.0,
            delr,
            270.0,
            delc,
        ]
        for icol in range(ncol_split)
    ]

    # set up 4 exchanges
    # tl-tr
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(lr_data),
        exgmnamea="tl",
        exgmnameb="tr",
        exchangedata=lr_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=useXT3D,
        filename="tltr.exg",
    )

    # bl-br
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(lr_data),
        exgmnamea="bl",
        exgmnameb="br",
        exchangedata=lr_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=useXT3D,
        filename="blbr.exg",
    )

    # tl-bl
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(tb_data),
        exgmnamea="tl",
        exgmnameb="bl",
        exchangedata=tb_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=useXT3D,
        filename="tlbl.exg",
    )

    # tr-br
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(tb_data),
        exgmnamea="tr",
        exgmnameb="br",
        exchangedata=tb_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=useXT3D,
        filename="trbr.exg",
    )

    return sim


def create_gwf_model(sim, mname, dis_params):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname, save_flows=True)

    dis_params = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=dis_params.nlay,
        nrow=dis_params.nrow,
        ncol=dis_params.ncol,
        delr=dis_params.delr,
        delc=dis_params.delc,
        xorigin=dis_params.xorigin,
        yorigin=dis_params.yorigin,
        top=dis_params.tops[0],
        botm=dis_params.tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=useXT3D,
        save_flows=True,
        icelltype=0,
        k=k11,
        k22=k22,
        angle1=k_angle,
    )

    # chd file
    left_chd = []
    right_chd = []
    if mname == "ref":
        left_chd = [[(0, irow, 0), h_left] for irow in range(nrow)]
        right_chd = [[(0, irow, ncol - 1), h_right] for irow in range(nrow)]
    elif mname == "tl" or mname == "bl":
        left_chd = [[(0, irow, 0), h_left] for irow in range(nrow_split)]
        right_chd = []
    elif mname == "tr" or mname == "br":
        left_chd = []
        right_chd = [[(0, irow, ncol_split - 1), h_right] for irow in range(nrow_split)]
    chd_data = left_chd + right_chd
    chd_spd = {0: chd_data}
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # well in top-left corner
    if mname == "ref" or mname == "tl":
        wel1 = flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=[[well_id, well_rate]],
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="WEL-1",
        )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname}.hds",
        budget_filerecord=f"{mname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )


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
    fpth = os.path.join(test.workspace, "ref.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()
    fpth = os.path.join(test.workspace, "ref.cbc")
    nlay, nrow, ncol = heads.shape
    qx, qy, qz = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, "tl.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_tl = hds.get_data()
    fpth = os.path.join(test.workspace, "tl.cbc")
    nlay, nrow, ncol = heads_tl.shape
    qx_tl, qy_tl, qz_tl = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, "tr.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_tr = hds.get_data()
    fpth = os.path.join(test.workspace, "tr.cbc")
    nlay, nrow, ncol = heads_tr.shape
    qx_tr, qy_tr, qz_tr = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, "bl.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_bl = hds.get_data()
    fpth = os.path.join(test.workspace, "bl.cbc")
    nlay, nrow, ncol = heads_bl.shape
    qx_bl, qy_bl, qz_bl = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, "br.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_br = hds.get_data()
    fpth = os.path.join(test.workspace, "br.cbc")
    nlay, nrow, ncol = heads_br.shape
    qx_br, qy_br, qz_br = qxqyqz(fpth, nlay, nrow, ncol)

    heads_top = np.append(heads_tl[0], heads_tr[0], axis=1)
    heads_bot = np.append(heads_bl[0], heads_br[0], axis=1)
    heads_merged = np.append(heads_top, heads_bot, axis=0)

    # compare heads
    maxdiff = np.amax(abs(heads - heads_merged))
    assert maxdiff < 10 * hclose_check, (
        f"Max. head diff. {maxdiff} should \
                     be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis-x
    qx_top = np.append(qx_tl[0], qx_tr[0], axis=1)
    qx_bot = np.append(qx_bl[0], qx_br[0], axis=1)
    qx_merged = np.append(qx_top, qx_bot, axis=0)

    maxdiff = np.amax(abs(qx - qx_merged))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (x) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis-y
    qy_top = np.append(qy_tl[0], qy_tr[0], axis=1)
    qy_bot = np.append(qy_bl[0], qy_br[0], axis=1)
    qy_merged = np.append(qy_top, qy_bot, axis=0)

    maxdiff = np.amax(abs(qy - qy_merged))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (y) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis-z
    qz_top = np.append(qz_tl[0], qz_tr[0], axis=1)
    qz_bot = np.append(qz_bl[0], qz_br[0], axis=1)
    qz_merged = np.append(qz_top, qz_bot, axis=0)

    maxdiff = np.amax(abs(qz - qz_merged))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (z) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # check budget error from .lst file
    for mname in ["ref", "tl", "tr", "bl", "br"]:
        fpth = os.path.join(test.workspace, f"{mname}.lst")
        for line in open(fpth):
            if line.lstrip().startswith("PERCENT"):
                cumul_balance_error = float(line.split()[3])
                assert abs(cumul_balance_error) < 0.00001, (
                    f"Cumulative balance error = {cumul_balance_error} for {mname}, "
                    "should equal 0.0"
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

"""
Test the exchange mover functionality to work both ways,
also in parallel. Use the following setup of two connected
DIS models with a stream flow crossing the boundary twice:


              left:                      right:
           .  .  .  .  .             .  .  .  .  .   1
 sfr_in => x  x  x  x  x             x  x  x  x  x   2
           .  .  .  .  .   gwfgwf    .  .  .  .  x   3
sfr_out <= x  x  x  x  x             x  x  x  x  x   4
           .  .  .  .  .             .  .  .  .  .   5
           1  2  3  4  5             6  7  8  9  10

The "single" model is also constructed as a reference.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwf_exgmvr02"]

nper = 1

lx = 10.0
ly = 5.0
lz = 1.0
nlay = 1
nrow = 5
ncol = 10
ncol_split = int(ncol / 2)
nper = 1
delc = ly / nrow
delr = lx / ncol
delz = lz / nlay
top = 0.0
botm = [top - (k + 1) * delz for k in range(nlay)]
Kh = 20.0
Kv = 20.0


def make_sfr_data(sfr_cells, ireach_offset=0):
    """generate package and connection data for a string of connected cells"""

    pak_data = []
    con_data = []

    rlen = delr
    rwid = delc
    rgrd = 1.0
    rtp = 0.0
    rbth = 0.1
    rhk = 0.01
    rman = 1.0
    ustrf = 1.0
    ndv = 0
    nc = len(sfr_cells)
    for ridx, cellid in enumerate(sfr_cells):
        irno = ridx + ireach_offset
        ncon = 2
        if ridx in [0, nc - 1]:
            ncon = 1
        t = (irno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, rman, ncon, ustrf, ndv)
        pak_data.append(t)

        if ridx == 0:  # first one only connected to the next
            c = (irno, -(irno + 1))
        elif ridx == nc - 1:  # last one only connected to the prev
            c = (irno, irno - 1)
        else:  # connect upstream and downstream
            c = (irno, irno - 1, -(irno + 1))
        con_data.append(c)

    return pak_data, con_data


def build_gwf(sim, gwf_type="single"):
    if gwf_type == "single":
        nc = ncol
    else:  # left or right
        nc = int(ncol / 2)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_type,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=nc,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=0,
        k=Kh,
        k33=Kv,
    )

    # add chd to right edge
    if gwf_type in ("single", "right"):
        chdlist = [[(0, irow, nc - 1), 0.0] for irow in range(nrow)]
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist,
            pname="chd_right",
        )

    left_sfr1 = [(0, 1, icol) for icol in range(ncol_split)]
    left_sfr2 = [(0, 3, icol - 1) for icol in range(ncol_split, 0, -1)]

    right_sfr = [(0, 1, icol) for icol in range(ncol_split)]
    right_sfr.append((0, 2, ncol_split - 1))
    right_sfr.extend([(0, 3, icol - 1) for icol in range(ncol_split, 0, -1)])

    # shift for single model:
    for t in right_sfr:
        print(t)
    right_sfr_single = [(t[0], t[1], t[2] + ncol_split) for t in right_sfr]

    package_data = []
    if gwf_type == "single":
        left_sfr1.extend(right_sfr_single)
        left_sfr1.extend(left_sfr2)
        package_data, conn_data = make_sfr_data(left_sfr1)
    elif gwf_type == "left":
        # these two are not connected
        package_data, conn_data = make_sfr_data(left_sfr1)
        package_data2, conn_data2 = make_sfr_data(
            left_sfr2, ireach_offset=len(package_data)
        )
        package_data.extend(package_data2)
        conn_data.extend(conn_data2)
    elif gwf_type == "right":
        package_data, conn_data = make_sfr_data(right_sfr)

    nreaches = len(package_data)

    if gwf_type in ("single", "left"):
        period_data = [
            (0, "INFLOW", 100.0),
        ]
    else:
        period_data = None

    if gwf_type != "single":
        mover = True
    else:
        mover = None

    sfr = flopy.mf6.modflow.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        mover=mover,
        stage_filerecord=f"{gwf_type}.sfr.stg",
        budget_filerecord=f"{gwf_type}.sfr.bud",
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=conn_data,
        perioddata=period_data,
        pname=f"sfr_{gwf_type}",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_type}.cbc",
        head_filerecord=f"{gwf_type}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "LAST")],
    )

    return gwf


def build_exchanges(sim):
    # add a gwf-gwf exchange
    gwfgwf_data = [
        (
            (0, irow, ncol_split - 1),
            (0, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            0.0,
            delr,
        )
        for irow in range(nrow)
    ]

    # GWF-GWF
    mvr_filerecord = "left-right.exg.mvr"
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea="left",
        exgmnameb="right",
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        dev_interfacemodel_on=False,
        filename="left-right.exg",
    )

    # simulation GWF-GWF Mover
    maxmvr, maxpackages = 2, 2
    mvrpack_sim = [["left", "sfr_left"], ["right", "sfr_right"]]
    mvrspd = [
        # connect left to right
        ["left", "sfr_left", ncol_split - 1, "right", "sfr_right", 0, "FACTOR", 1.00],
        # connect right to left
        [
            "right",
            "sfr_right",
            2 * ncol_split,
            "left",
            "sfr_left",
            ncol_split,
            "FACTOR",
            1.00,
        ],
    ]

    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=maxmvr,
        print_flows=True,
        maxpackages=maxpackages,
        packages=mvrpack_sim,
        perioddata=mvrspd,
        filename=mvr_filerecord,
    )


def build_simulation(idx, sim_ws, sim_type="single"):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=sim_ws)

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper)

    # Flow solver
    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="simple",
        print_option="ALL",
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
    )

    if sim_type == "single":
        gwf_types = ("single",)
    else:
        gwf_types = ("left", "right")
    for gwf_type in gwf_types:
        gwf = build_gwf(sim, gwf_type=gwf_type)

    if sim_type != "single":
        build_exchanges(sim)

    return sim


def build_models(idx, test):
    sim_ws = test.workspace / "mf6"
    sim_base = build_simulation(idx, sim_ws)
    sim = build_simulation(idx, test.workspace, sim_type="split")
    return sim, sim_base


def check_output(idx, test):
    # base simulations stage
    ws = test.workspace
    fpth = ws / "mf6/single.sfr.stg"
    single_stage_obj = flopy.utils.HeadFile(fpth, text="STAGE")
    single_stage = single_stage_obj.get_data().squeeze()

    fpth = ws / "left.sfr.stg"
    stage_obj = flopy.utils.HeadFile(fpth, text="STAGE")
    v = stage_obj.get_data().squeeze()
    assert np.allclose(single_stage[0:ncol_split], v[0:ncol_split]), (
        "sfr left (segment I) stages are not equal"
    )
    assert np.allclose(single_stage[3 * ncol_split + 1 :], v[ncol_split:]), (
        "sfr left (segment II) stages are not equal"
    )

    fpth = ws / "right.sfr.stg"
    stage_obj = flopy.utils.HeadFile(fpth, text="STAGE")
    v = stage_obj.get_data().squeeze()
    assert np.allclose(single_stage[ncol_split : 3 * ncol_split + 1], v), (
        "sfr right stages are not equal"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()

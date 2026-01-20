"""
test to evaluate Newton-Raphson solution for a single column transient
dry multi-aquifer well problem. Developed to address issue #546
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ("maw_09a", "maw_09b", "maw_09c", "maw_09d")
dis_option = ("dis", "dis", "disv", "disv")
flow_correction = (None, True, None, True)

nlay = 3
nrow = 1
ncol = 1

vertices = [
    (0, 0, 0),
    (1, 0, 1),
    (2, 1, 1),
    (3, 1, 0),
]
cell2d = [
    (0, 0.5, 0.5, 5, 0, 1, 2, 3, 0),
]

delc = 1.0
delr = 1.0
gwfarea = delr * delc

top = 30.0
botm = [20, 10, 0]
zelevs = [top] + botm
maw_bot = 15.0

strt = 4.5
maw_strt = 30.0

Kh = 10
Kv = 10.0
radius = np.sqrt(1.0 / np.pi)


def build_models(idx, test):
    dvclose, rclose, relax = 1e-9, 1e-9, 1.0

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="summary",
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[(20.0, 50, 1.1)])

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="complex",
        outer_dvclose=dvclose,
        inner_dvclose=dvclose,
        rcloserecord=f"{rclose} strict",
        relaxation_factor=relax,
    )

    # create gwf model
    gwfname = "gwf_" + name

    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions=newtonoptions, print_input=True
    )

    if dis_option[idx] == "dis":
        dis = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
        cellids = [
            (0, 0, 0),
            (1, 0, 0),
        ]
        gwf_obs = [
            ["C1", "HEAD", (0, 0, 0)],
            ["C2", "HEAD", (1, 0, 0)],
            ["C3", "HEAD", (2, 0, 0)],
        ]
    else:
        dis = flopy.mf6.ModflowGwfdisv(
            gwf,
            nlay=nlay,
            ncpl=nrow,
            nvert=4,
            vertices=vertices,
            cell2d=cell2d,
            top=top,
            botm=botm,
        )
        cellids = [
            (0, 0),
            (1, 0),
        ]
        gwf_obs = [
            ["C1", "HEAD", (0, 0)],
            ["C2", "HEAD", (1, 0)],
            ["C3", "HEAD", (2, 0)],
        ]

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=Kh,
        k33=Kv,
    )
    # gwf observations
    opth = f"{gwfname}.gwf.obs"
    obsdata = {f"{gwfname}.gwf.obs.csv": gwf_obs}
    flopy.mf6.ModflowUtlobs(
        gwf, filename=opth, digits=20, print_input=True, continuous=obsdata
    )

    # storage
    sto = flopy.mf6.ModflowGwfsto(gwf, ss=0.0, sy=1.0, transient=True, iconvert=1)

    # <ifno> <radius> <bottom> <strt> <condeqn> <ngwfnodes>
    mawpackagedata = flopy.mf6.ModflowGwfmaw.packagedata.empty(gwf, maxbound=1)
    mawpackagedata["radius"] = radius
    mawpackagedata["bottom"] = maw_bot
    mawpackagedata["strt"] = maw_strt
    mawpackagedata["condeqn"] = "specified"
    mawpackagedata["ngwfnodes"] = 2

    # <ifno> <icon> <cellid(ncelldim)> <scrn_top> <scrn_bot> <hk_skin> <radius_skin>
    mawconnectiondata = flopy.mf6.ModflowGwfmaw.connectiondata.empty(gwf, maxbound=2)
    mawconnectiondata["icon"] = [0, 1]
    mawconnectiondata["cellid"] = cellids
    mawconnectiondata["scrn_top"] = 100.0
    mawconnectiondata["scrn_bot"] = 0.0
    mawconnectiondata["hk_skin"] = 1.0
    mawconnectiondata["radius_skin"] = -999

    mbin = f"{gwfname}.maw.bin"
    mbud = f"{gwfname}.maw.bud"
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        flow_correction=flow_correction[idx],
        head_filerecord=mbin,
        budget_filerecord=mbud,
        packagedata=mawpackagedata,
        connectiondata=mawconnectiondata,
    )
    opth = f"{gwfname}.maw.obs"
    obsdata = {
        f"{gwfname}.maw.obs.csv": [
            ("whead", "head", (0,)),
        ]
    }
    maw.obs.initialize(filename=opth, digits=20, print_input=True, continuous=obsdata)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
        printrecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwfname = "gwf_" + name

    # get well observations
    fname = f"{gwfname}.gwf.obs.csv"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    gobs = np.genfromtxt(fname, delimiter=",", names=True)

    # get well observations
    fname = f"{gwfname}.maw.obs.csv"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    wobs = np.genfromtxt(fname, delimiter=",", names=True)["WHEAD"]

    # calculate volume of water and make sure it is conserved
    # volume comparisons can be made based on saturated thickness because
    # the cell area and well area are both equal to 1
    v0 = (maw_strt - 10.0) + (strt - 0.0)
    for i, w in enumerate(wobs):
        vg = 0.0
        for jdx, tag in enumerate(("C1", "C2", "C3")):
            g = gobs[tag][i]
            ctop = zelevs[jdx]
            cbot = zelevs[jdx + 1]
            if g > ctop:
                d = ctop - cbot
            elif g > cbot:
                d = g - cbot
            else:
                d = 0.0
            vg += d
        vw = w - 10.0
        vt = vw + vg

        # write a message
        msg = f"{i}\n  well volume: {vw} "
        msg += f"\n  groundwater volume: {vg}"
        msg += f"\n  total volume: {vt}"
        print(msg)

        # evaluate results
        msg = f"total volume {vt} not equal to initial volume {v0}"
        msg += f"\nwell volume: {vw}\ngroundwater volume: {vg}"
        assert np.allclose(vt, v0), msg

    # Check final well head
    well_head = wobs[-1]
    assert np.allclose(well_head, 17.25), (
        f"final simulated maw head ({well_head}) does not equal 17.25."
    )

    fname = gwfname + ".hds"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    hobj = flopy.utils.HeadFile(fname)
    head = hobj.get_alldata()[1:]

    # compare the maw-gwf flows with the gwf-maw flows
    fname = gwfname + ".maw.bud"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    mbud = flopy.utils.CellBudgetFile(fname, precision="double")
    maw_gwf = mbud.get_data(text="GWF")

    fname = gwfname + ".cbc"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)
    gbud = flopy.utils.CellBudgetFile(fname, precision="double")
    gwf_maw = gbud.get_data(text="MAW")

    assert len(maw_gwf) == len(gwf_maw), "number of budget records not equal"

    for istp, (ra_maw, ra_gwf) in enumerate(zip(maw_gwf, gwf_maw)):
        for i in range(ra_maw.shape[0]):
            qmaw = ra_maw[i]["q"]
            qgwf = ra_gwf[i]["q"]
            msg = f"step {istp} record {i} comparing qmaw with qgwf: {qmaw} {qgwf}"
            print(msg)
            assert np.allclose(qmaw, -qgwf), msg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        compare="mf6_regression",
    )
    test.run()

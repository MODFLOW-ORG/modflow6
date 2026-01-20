"""
Based on sft01 gwf model, but split into two gwf models test gwf-gwf and
mvr. The single model is run as the regression model

The final split model look like:

      flow1                        flow2
 sfr  1 2 3 4 5 6 7  gwfgwf-mvr => 1 2 3 4 5 6 7
      -------------                -------------
 gwf  1 2 3 4 5 6 7  gwfgwf     => 1 2 3 4 5 6 7
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = ["gwf_exgmvr01"]

# properties for single model combination
lx = 14.0
lz = 1.0
nlay = 1
nrow = 1
ncol = 14
nper = 1
delc = 1.0
delr = lx / ncol
delz = lz / nlay
top = 0.0
botm = [top - (k + 1) * delz for k in range(nlay)]
Kh = 20.0
Kv = 20.0


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


def build_gwf(sim, gwf_type="single"):
    if gwf_type == "single":
        nc = ncol
    else:
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
        chdlist = [
            [(0, 0, nc - 1), 0.0],
        ]
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=chdlist,
            pname="chd_right",
        )

    # inject water into left edge
    if gwf_type in ("single", "left"):
        wellist = [
            [(0, 0, 0), 1.0],
        ]
        wel = flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data=wellist,
            pname="well_left",
        )

    # pak_data = [<rno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #             <man> <ncon> <ustrf> <ndv> [<aux(naux)>] [<boundname>]]
    rlen = delr
    rwid = delc
    rgrd = 1.0
    rtp = 0.0
    rbth = 0.1
    rhk = 0.01
    rman = 1.0
    ustrf = 1.0
    ndv = 0
    pak_data = []
    for irno in range(nc):
        ncon = 2
        if irno in [0, nc - 1]:
            ncon = 1
        cellid = (0, 0, irno)
        t = (irno, cellid, rlen, rwid, rgrd, rtp, rbth, rhk, rman, ncon, ustrf, ndv)
        pak_data.append(t)

    con_data = []
    for irno in range(nc):
        if irno == 0:
            t = (irno, -(irno + 1))
        elif irno == nc - 1:
            t = (irno, irno - 1)
        else:
            t = (irno, irno - 1, -(irno + 1))
        con_data.append(t)

    if gwf_type in ("single", "left"):
        p_data = [
            (0, "INFLOW", 1.0),
        ]
    else:
        p_data = None

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
        nreaches=nc,
        packagedata=pak_data,
        connectiondata=con_data,
        perioddata=p_data,
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
            (0, 0, int(ncol / 2) - 1),
            (0, 0, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            0.0,
            delr,
        )
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
    maxmvr, maxpackages = 1, 2
    mvrpack_sim = [["left", "sfr_left"], ["right", "sfr_right"]]
    mvrspd = [
        ["left", "sfr_left", int(ncol / 2) - 1, "right", "sfr_right", 0, "FACTOR", 1.00]
    ]

    gwfgwf.mvr.initialize(
        modelnames=True,
        maxmvr=maxmvr,
        print_flows=True,
        budgetcsv_filerecord="left-right.mvr.bud.csv",
        maxpackages=maxpackages,
        packages=mvrpack_sim,
        perioddata=mvrspd,
        filename=mvr_filerecord,
    )


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

    stage = single_stage.copy()

    i1 = int(ncol / 2)
    fpth = ws / "left.sfr.stg"
    stage_obj = flopy.utils.HeadFile(fpth, text="STAGE")
    v = stage_obj.get_data().squeeze()
    stage[:i1] = v[:]

    fpth = ws / "right.sfr.stg"
    stage_obj = flopy.utils.HeadFile(fpth, text="STAGE")
    v = stage_obj.get_data().squeeze()
    stage[i1:] = v[:]

    assert np.allclose(single_stage, stage), "sfr stages are not equal"

    # check mover CSV output
    fpth = ws / "left-right.mvr.bud.csv"
    df = pd.read_csv(fpth)

    for diff in df["PERCENT_DIFFERENCE"].iloc:
        assert abs(diff) < 1e-05, "no mass balance on the exchange mover"


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

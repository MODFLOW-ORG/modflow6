# Test evap in SFR reaches (no interaction with gwf)

import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

HDRY, HNOFLO = -1e30, 1e30

cases = ["sfr-inactive02"]


def get_model(ws, name, array_input=False):
    # Base simulation and model name and workspace

    length_units = "meters"
    time_units = "seconds"

    nrow = 1
    ncol = 1
    nlay = 1
    delr = delc = 1.0

    nper = 2
    tds_spd = [
        (1.0, 1, 1.0),
        (1.0, 1, 1.0),
    ]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units=time_units,
        nper=nper,
        perioddata=tds_spd,
    )
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=-100.0,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,  # >0 means saturated thickness varies with computed head
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    if array_input:
        # if False:
        bhead = np.full(nlay * nrow * ncol, DNODATA, dtype=float)
        cond = np.full(nlay * nrow * ncol, DNODATA, dtype=float)
        bhead[0] = 1.0
        cond[0] = 1e6
        flopy.mf6.ModflowGwfghbg(gwf, maxbound=1, bhead=bhead, cond=cond)
    else:
        flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[((0, 0, 0), 1.0, 1e6)])

    # sfr data
    nreaches = 4

    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #        <man> <ncon> <ustrf> <ndv>
    package_data = [
        (0, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, 1.0, 0.001, 1, 0.0, 0),
        (1, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, 1.0, 0.001, 1, 1.0, 0),
        (2, (-1, -1, -1), delr, 1.0, 1e-3, 0.0, 1.0, 1.0, 0.001, 0, 0.0, 0),
        (3, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, 1.0, 0.001, 0, 0.0, 0),
    ]
    connection_data = [
        (0, -1),
        (1, 0),
        (2,),
        (3,),
    ]

    sfr_spd = {
        0: [
            (0, "inflow", 1.0),
            (0, "rainfall", 1.0),
            (0, "evaporation", 1.0),
            (0, "runoff", 1.0),
            (1, "rainfall", 1.0),
            (1, "evaporation", 1.0),
            (1, "runoff", 1.0),
        ],
        1: [
            (1, "status", "inactive"),
        ],
    }

    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwf", "sfr", (0,)),
            ("outflow", "ext-outflow", (0,)),
            ("depth", "depth", (0,)),
        ],
        "filename": name + ".sfr.obs",
    }

    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        stage_filerecord=f"{name}.sfr.stg",
        budget_filerecord=f"{name}.sfr.cbc",
        length_conversion=1.0,
        time_conversion=1.0,
        mover=True,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=sfr_spd,
        observations=sfr_obs,
        pname="SFR-1",
    )

    flopy.mf6.ModflowGwfmvr(
        gwf,
        print_input=True,
        print_flows=True,
        maxmvr=1,
        maxpackages=1,
        packages=["SFR-1"],
        perioddata={
            1: [
                ("SFR-1", 0, "SFR-1", 3, "FACTOR", 1.0),
            ],
        },
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("head", "all"), ("budget", "all")],
    )

    return sim


def check_output(ws, name):
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwf = sim.get_model()
    sfr = gwf.get_package("SFR-1")
    stage = sfr.output.stage().get_alldata().squeeze()
    for idx in (0, 1):
        assert stage[idx, 2] == HDRY, (
            f"reach 3 stage is not HDRY in stress period {idx + 1}"
        )
    assert stage[1, 1] == HNOFLO, "reach 4 stage is not HNOFLO in stress period 2"

    bobj = sfr.output.budget()
    data_names = (
        "GWF",
        "RAINFALL",
        "EVAPORATION",
        "RUNOFF",
        "EXT-INFLOW",
        "EXT-OUTFLOW",
        "STORAGE",
        "FROM-MVR",
        "TO-MVR",
    )
    for name in data_names:
        v = bobj.get_data(text=name, totim=2.0)[0]
        assert v["q"][1] == 0.0, f"{name} flow for reach 2 is not zero ({v['q'][1]})"

    # skip GWF for reach 3 since it is not connected
    # to a GWF cell (data_names[0])
    for totim in (1.0, 2.0):
        for name in data_names[1:]:
            v = bobj.get_data(text=name, totim=totim)[0]
            assert v["q"][2] == 0.0, (
                f"{name} flow for reach 3 is not zero ({v['q'][2]})"
            )

    v = bobj.get_data(text="FLOW-JA-FACE", totim=1.0)[0]
    node = v["node"]
    node2 = v["node2"]
    conn = (node[0], node2[0])
    conn2 = (node2[1], node[1])
    assert conn == conn2, "FLOW-JA-FACE connectivity is not symmetric"

    q = v["q"]
    assert np.allclose(q[0], -q[1]), (
        f"FLOW-JA-FACE for 1-2 ({q[0]}) not equal to 2-1 ({-q[1]})"
    )

    fa = v["FLOW-AREA"]
    assert np.allclose(fa[0], fa[1]), (
        f"FLOW-AREA for 1-2 ({fa[0]}) not equal to 2-1 ({fa[1]})"
    )

    v = bobj.get_data(text="TO-MVR", totim=2.0)[0]
    v2 = bobj.get_data(text="FROM-MVR", totim=2.0)[0]
    assert np.allclose(v["q"][0], -v2["q"][3]), (
        f"TO-MVR for reach 1 ({v['q'][0]}) not equal "
        + f"to FROM-MVR for reach 4 ({v2['q'][3]})"
    )


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name)

    # build comparison array_input model
    ws = os.path.join(test.workspace, "mf6")
    mc = get_model(ws, name, array_input=True)

    return sim, mc


def check_outputs(idx, test):
    name = cases[idx]

    # check output MODFLOW 6 files
    ws = test.workspace
    check_output(ws, name)

    # check output comparison array_input model
    ws = os.path.join(test.workspace, "mf6")
    check_output(ws, name)


@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_outputs(idx, t),
        compare="mf6",
    )
    test.run()

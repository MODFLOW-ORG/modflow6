"""
Test for lake package outlet observations.

The test evaluates the total outlet flow calculated using lake boundname
with EXT-OUTFLOW and OUTLET flow for all three outlets.
"""

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["lakoutlet_obs"]


def build_models(idx, test):
    nlay = 2
    nrow = 3
    ncol = 3
    nper = 1
    delc = 1000.0
    delr = 1000.0
    top = 5.0
    botm = [-5.0, -10.0]

    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]

    Kh = 1.0

    tdis_rc = [(1.0, 1, 1.0)]

    name = cases[idx]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)  # , newtonoptions="newton")

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="simple",
        inner_dvclose=1e-9,
        outer_dvclose=1e-8,
        # linear_acceleration="bicgstab",
    )

    # number of columns to be a lake for layer 1, 2, , ... len(lakend)
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[0, 1, 1] = 0
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
    strt = 0.0
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=Kh,
    )

    sy = 0.3
    ss = 0.0
    sto = flopy.mf6.ModflowGwfsto(gwf, sy=sy, ss=ss, iconvert=1)

    lake_connect = [(0, 1), (1, 0), (1, 2), (2, 1)]
    nlakeconn = len(lake_connect) + 1

    # pak_data = [ifno, strt, nlakeconn]
    pak_data = [(0, strt, nlakeconn, "LAKE1")]

    bedleak = DNODATA
    belev = botm[0]
    con_data = [
        (0, idx, (0, i, j), "HORIZONTAL", bedleak, belev, top, delr / 2.0, delr)
        for idx, (i, j) in enumerate(lake_connect)
    ]
    con_data.append((0, 4, (1, 1, 1), "VERTICAL", bedleak, belev, top, 0.0, 0.0))

    # outlet data
    outlet_data = [
        (0, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0),
        (1, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0),
        (2, 0, -1, "specified", 0.0, 0.0, 0.0, 0.0),
    ]
    noutlets = len(outlet_data)

    # period data
    p_data = [
        (0, "RATE", -3.0),
        (1, "RATE", -2.0),
        (2, "RATE", -1.0),
        # (0, "INFLOW", 6.0),
    ]

    lak_obs = {
        f"{name}.lak.obs.csv": [
            ("total", "ext-outflow", "LAKE1"),
            ("o1", "outlet", (0,)),
            ("o2", "outlet", (1,)),
            ("o3", "outlet", (2,)),
        ],
    }

    lak = flopy.mf6.modflow.ModflowGwflak(
        gwf,
        surfdep=0.0,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_stage=True,
        boundnames=True,
        nlakes=len(pak_data),
        ntables=0,
        noutlets=noutlets,
        outlets=outlet_data,
        packagedata=pak_data,
        pname="LAK-1",
        connectiondata=con_data,
        perioddata=p_data,
        observations=lak_obs,
    )

    # idx = np.where(idomain == 1)
    # chdspd = []
    # for k, i, j in zip(*idx):
    #     chdspd.append(((k, i, j), 0.0))
    # chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=chdspd)

    rech = 6.0 / (8.0 * delr * delc)
    rch = flopy.mf6.modflow.ModflowGwfrcha(
        gwf,
        print_flows=True,
        save_flows=True,
        recharge=rech,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    obs_values = flopy.utils.Mf6Obs(test.workspace / f"{name}.lak.obs.csv")
    test_values = {
        "TOTAL": [-6.0],
        "O1": [-3.0],
        "O2": [-2.0],
        "O3": [-1.0],
    }
    for key, value in test_values.items():
        assert np.array_equal(obs_values.get_data(obsname=key)[key], value), (
            f"failed comparison for '{key}' observation"
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

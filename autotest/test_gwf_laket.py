"""Test for checking lak evaporation."""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "gwf_laket01",
    "gwf_laket02",
    "gwf_laket03",
]

bedleak = [0.0, 1.0, 1.0]
strt = [5.0, 6.0, 4.0]
lakestage = [6.0, 5.0, 6.0]


def get_model(idx, ws):
    name = cases[idx]
    nlay = 1
    nrow = 1
    ncol = 1
    nper = 2
    delc = 1.0
    delr = 1.0
    delz = 10.0
    top = 5.0
    botm = [top - (k + 1) * delz for k in range(nlay)]

    perlen = [11.0, 11.0]
    nstp = [11, 11]
    tsmult = [1.0, 1.0]

    Kh = 1.0
    Kv = 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        complexity="simple",
        linear_acceleration="BICGSTAB",
    )

    # create gwf model
    gwfname = name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="NEWTON")

    # discretization package
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

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt[idx])

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_specific_discharge=False,
        icelltype=1,
        k=Kh,
        k33=Kv,
    )

    sy = 0.3
    ss = 1e-5
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        ss_confined_only=False,
        sy=sy,
        ss=ss,
        iconvert=1,
    )

    lake_vconnect = [
        (0, 0, 0),
    ]
    nlakeconn = 1

    # pak_data = [ifno, strt, nlakeconn]
    pak_data = [(0, lakestage[idx], nlakeconn)]

    belev = top
    con_data = [
        (0, i, cellid, "VERTICAL", bedleak[idx], belev, -99, -99, -99)
        for i, cellid in enumerate(lake_vconnect)
    ]

    # period data
    p_data = {
        0: [
            (0, "EVAPORATION", 0.1),
        ],
        1: [
            (0, "RAINFALL", 0.2),
        ],
    }

    # note: for specifying lake number, use fortran indexing!
    fname = f"{gwfname}.lak.obs.csv"
    lak_obs = {
        fname: [
            ("lakestage", "stage", (0,)),
            ("lakevolume", "volume", (0,)),
            ("evap", "evaporation", (0,)),
        ],
    }

    lak = flopy.mf6.modflow.ModflowGwflak(
        gwf,
        surfdep=0.0,
        print_input=True,
        print_flows=True,
        print_stage=True,
        nlakes=len(pak_data),
        ntables=0,
        packagedata=pak_data,
        pname="LAK-1",
        connectiondata=con_data,
        perioddata=p_data,
        observations=lak_obs,
        package_convergence_filerecord=f"{gwfname}.lak.convergence.csv",
    )

    if idx > 0:
        ghb = flopy.mf6.modflow.ModflowGwfghb(
            gwf,
            stress_period_data=[
                (0, 0, 0, strt[idx], 1000.0),
            ],
        )
        # chd = flopy.mf6.modflow.ModflowGwfchd(
        #     gwf,
        #     stress_period_data=[
        #         (0, 0, 0, strt[idx]),
        #     ],
        # )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("BUDGET", "ALL"), ("HEAD", "ALL")],
    )

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, f"{test.name}.lak.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    dtype = [("time", float), ("stage", float), ("evap", float)]
    obs = {
        0: np.array(
            [
                (1.000000000000, 5.9000000000000004, -0.1),
                (2.000000000000, 5.8000000000000007, -0.1),
                (3.000000000000, 5.7000000000000011, -0.1),
                (4.000000000000, 5.6000000000000014, -0.1),
                (5.000000000000, 5.5000000000000018, -0.1),
                (6.000000000000, 5.4000000000000021, -0.1),
                (7.000000000000, 5.3000000000000025, -0.1),
                (8.000000000000, 5.2000000000000028, -0.1),
                (9.000000000000, 5.1000000000000032, -0.1),
                (10.00000000000, 5.0000000000000036, -0.1),
                (11.00000000000, 5.0000000000000000, 0.0),
                (12.00000000000, 5.0999999999999996, -0.1),
                (13.00000000000, 5.1999999999999993, -0.1),
                (14.00000000000, 5.2999999999999989, -0.1),
                (15.00000000000, 5.3999999999999986, -0.1),
                (16.00000000000, 5.4999999999999982, -0.1),
                (17.00000000000, 5.5999999999999979, -0.1),
                (18.00000000000, 5.6999999999999975, -0.1),
                (19.00000000000, 5.7999999999999972, -0.1),
                (20.00000000000, 5.8999999999999968, -0.1),
                (21.00000000000, 5.9999999999999964, -0.1),
                (22.00000000000, 6.0999999999999961, -0.1),
            ],
            dtype=dtype,
        ),
        1: np.array(
            [
                (1.000000000000, 5.0571204119063822, -0.1),
                (2.000000000000, 5.1060819299607383, -0.1),
                (3.000000000000, 5.1480499445118104, -0.1),
                (4.000000000000, 5.1840233847838952, -0.1),
                (5.000000000000, 5.2148584962098710, -0.1),
                (6.000000000000, 5.2412892209184587, -0.1),
                (7.000000000000, 5.2639446671391816, -0.1),
                (8.000000000000, 5.2833640833345035, -0.1),
                (9.000000000000, 5.3000096934757943, -0.1),
                (10.00000000000, 5.3142776989703462, -0.1),
                (11.00000000000, 5.3265077091090305, -0.1),
                (12.00000000000, 5.5084234768802327, -0.1),
                (13.00000000000, 5.6643549898781638, -0.1),
                (14.00000000000, 5.7980137542817234, -0.1),
                (15.00000000000, 5.9125811368359047, -0.1),
                (16.00000000000, 6.0107840882389674, -0.1),
                (17.00000000000, 6.0949600504480319, -0.1),
                (18.00000000000, 6.1671125928369053, -0.1),
                (19.00000000000, 6.2289591014666277, -0.1),
                (20.00000000000, 6.2819716565762063, -0.1),
                (21.00000000000, 6.3274120712659299, -0.1),
                (22.00000000000, 6.3663619253694703, -0.1),
            ],
            dtype=dtype,
        ),
        2: np.array(
            [
                (1.000000000000, 5.7714285714285714, -0.1),
                (2.000000000000, 5.5755102040816329, -0.1),
                (3.000000000000, 5.4075801749271140, -0.1),
                (4.000000000000, 5.2636401499375260, -0.1),
                (5.000000000000, 5.1402629856607369, -0.1),
                (6.000000000000, 5.0345111305663464, -0.1),
                (7.000000000000, 5.0000000000000000, -0.0345111305663464),
                (8.000000000000, 5.0000000000000000, 0.0),
                (9.000000000000, 5.0000000000000000, 0.0),
                (10.00000000000, 5.0000000000000000, 0.0),
                (11.00000000000, 5.0000000000000000, 0.0),
                (12.00000000000, 5.0857142857142854, -0.1),
                (13.00000000000, 5.1591836734693874, -0.1),
                (14.00000000000, 5.2221574344023320, -0.1),
                (15.00000000000, 5.2761349437734273, -0.1),
                (16.00000000000, 5.3224013803772232, -0.1),
                (17.00000000000, 5.3620583260376202, -0.1),
                (18.00000000000, 5.3960499937465318, -0.1),
                (19.00000000000, 5.4251857089255990, -0.1),
                (20.00000000000, 5.4501591790790851, -0.1),
                (21.00000000000, 5.4715650106392157, -0.1),
                (22.00000000000, 5.4899128662621850, -0.1),
            ],
            dtype=dtype,
        ),
    }

    if idx in (0, 1, 2):
        evap_compare = np.allclose(obs[idx]["evap"], tc["EVAP"])
        stage_compare = np.allclose(obs[idx]["stage"], tc["LAKESTAGE"])
    else:
        evap_compare = True
        stage_compare = True

    test.success = True
    msg = ""
    if not evap_compare:
        test.success = False
        msg += " Lake evaporation comparison failed."
    if not stage_compare:
        test.success = False
        msg += " Lake stage comparison failed."
    assert test.success, msg


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

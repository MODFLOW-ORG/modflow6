import os

import flopy
import numpy as np
import pytest
from cross_section_functions import calculate_rectchan_mannings_discharge
from framework import TestFramework

paktest = "sfr"
cases = [
    "sfr_npt03a",
    "sfr_npt03b",
    "sfr_npt03c",
    "sfr_npt03d",
    "sfr_npt03e",
    "sfr_npt03f",
    "sfr_npt03g",
]

simulated_depths = (
    0.5,
    1.0,
    1.1,
    2.0,
    3.0,
    10.0,
    100.0,
)

# temporal discretization
nper = 1
tdis_rc = []
for n in range(nper):
    tdis_rc.append((1.0, 1, 1.0))

# spatial discretization data
nlay, nrow, ncol = 1, 1, 1
delr, delc = 100.0, 100.0
top = 0.0
botm = -10.0
strt = 0.0


# static sfr data
rlen = 50.0
conversion_fact = 1.0
roughness = 0.01
rbth = 1.0
rhk = 0.0
slope = 0.001
ustrf = 1.0
ndv = 0
inflow = 1000.0

np_data = {
    "x": np.array([0.0, (1.0 / 3.0), (2.0 / 3.0), 1.0], dtype=float),
    "h": np.array([0.0, 0.0, 0.0, 0.0], dtype=float),
    "r": np.array([10.0, 1.0, 10.0, 10.0], dtype=float),
}

#  Cross section depiction
#
#    |           |           |           |
#    |  (left)   | (channel) |  (right)  |
#    |   10.0    |    1.0    |   10.0    | <- "MANFRACTION"
#    |           |           |           |
#    +-----------+-----------+-----------+  y: 0, 0, 0, 0
# x: 0          1/3         2/3          1
#


def build_model(idx, ws, base=False):
    if base:
        ws = os.path.join(ws, "mf6")

    # build MODFLOW 6 files
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=False,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf)

    # chd files
    # chd data
    spd = [
        [(0, 0, 0), 0.0],
    ]
    chd = flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd, pname="chd-1")

    # sfr data
    if base:
        nreaches = 6
        rwid = 10.0
        connection_data = [(0, -1), (1, 0), (2, -3), (3, 2), (4, -5), (5, 4)]
    else:
        nreaches = 2
        rwid = 30.0
        connection_data = [
            (0, -1),
            (1, 0),
        ]
    nconn = 1
    package_data = []
    for irch in range(nreaches):
        rfact = 1.0
        if base:
            if irch in (0, 1, 4, 5):
                rfact = 10.0
            if irch in (0, 2, 4):
                obsname = "upstream"
            else:
                obsname = "downstream"
        else:
            if irch == 0:
                obsname = "upstream"
            else:
                obsname = "downstream"
        rp = [
            irch,
            "none",
            rlen,
            rwid,
            slope,
            top,
            rbth,
            rhk,
            roughness * rfact,
            nconn,
            ustrf,
            ndv,
            obsname,
        ]
        package_data.append(rp)

    depth = simulated_depths[idx]
    qleft = calculate_rectchan_mannings_discharge(
        conversion_fact, roughness * 10.0, slope, 10.0, depth
    )
    qchannel = calculate_rectchan_mannings_discharge(
        conversion_fact, roughness, slope, 10.0, depth
    )
    qright = calculate_rectchan_mannings_discharge(
        conversion_fact, roughness * 10.0, slope, 10.0, depth
    )
    qtotal = qleft + qchannel + qright

    # create period data, observationsm and cross-section files
    fname = f"{name}.sfr.obs"
    if base:
        period_data = {
            0: [
                (0, "inflow", qleft),
                (2, "inflow", qchannel),
                (4, "inflow", qright),
            ]
        }
        sfr_obs = {
            f"{fname}.csv": [
                ("outflow_upstream", "outflow", "upstream"),
                ("outflow_downstream", "ext-outflow", "downstream"),
                ("depth_left", "depth", (1,)),
                ("depth_channel", "depth", (3,)),
                ("depth_right", "depth", (5,)),
                ("depth_upstream", "depth", (2,)),
            ]
        }
    else:
        sfr_tab = f"{name}.{n:02d}.sfr.tab"
        pname = f"sfrtab{n:02d}"
        stations = np_data["x"]
        depths = np_data["h"]
        roughs = np_data["r"]
        table = [[x, d, r] for x, d, r in zip(stations, depths, roughs)]
        t = flopy.mf6.ModflowUtlsfrtab(
            gwf,
            nrow=stations.shape[0],
            ncol=3,
            table=table,
            filename=sfr_tab,
            pname=pname,
        )
        t.write()
        period_data = {
            0: [
                (0, "inflow", qtotal),
                (0, "cross_section", sfr_tab),
                (1, "cross_section", sfr_tab),
            ]
        }
        sfr_obs = {
            f"{fname}.csv": [
                ("outflow_upstream", "outflow", "upstream"),
                ("outflow_downstream", "ext-outflow", "downstream"),
                ("depth", "depth", (1,)),
                ("depth_upstream", "depth", (0,)),
            ]
        }

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        boundnames=True,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=period_data,
        pname="sfr-1",
    )
    sfr.obs.initialize(
        filename=fname,
        digits=25,
        print_input=True,
        continuous=sfr_obs,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[
            ("BUDGET", "ALL"),
        ],
    )

    return sim


def build_models(idx, test):
    sim = build_model(idx, test.workspace)
    mc = build_model(idx, test.workspace, base=True)
    return sim, mc


def check_output(idx, test):
    obs_pth0 = os.path.join(test.workspace, f"{test.name}.sfr.obs.csv")
    obs0 = np.genfromtxt(obs_pth0, names=True, delimiter=",")

    obs_pth1 = os.path.join(test.workspace, "mf6", f"{test.name}.sfr.obs.csv")
    obs1 = np.genfromtxt(obs_pth1, names=True, delimiter=",")

    q0 = obs0["OUTFLOW_DOWNSTREAM"]
    q1 = obs1["OUTFLOW_DOWNSTREAM"]
    assert np.allclose(q0, q1), f"downstream outflows not equal ('{test.name}')"

    d0 = obs0["DEPTH_UPSTREAM"]
    d1 = obs1["DEPTH_UPSTREAM"]
    assert np.allclose(d0, d1), f"upstream depths are not equal ('{test.name}')"


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

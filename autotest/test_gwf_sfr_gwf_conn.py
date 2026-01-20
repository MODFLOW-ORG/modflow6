import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = [
    "sfr_dis",
    "sfr_dis_fail",
    "sfr_dis_none",
    "sfr_disv",
    "sfr_disv_fail",
    "sfr_disv_none",
    "sfr_disu",
    "sfr_disu_fail",
    "sfr_disu_none",
]
dis_types = [
    "dis",
    "dis",
    "dis",
    "disv",
    "disv",
    "disv",
    "disu",
    "disu",
    "disu",
]

# spatial discretization data
nlay, nrow, ncol = 1, 1, 1
delr, delc = 100.0, 100.0
top = 0.0
botm = -10.0
strt = 0.0

# spatial discretization data for disv and disu
vertices = [(0, 0.0, 0.0), (1, 0.0, delc), (2, delr, delc), (3, delr, 0.0)]
cell2d = [(0, delr / 2.0, delc / 2.0, 4, 0, 1, 2, 3)]

# sfr data
nreaches = 10
rlen = 10.0
rwid = 10.0
roughness = 0.001
rbth = 1.0
rhk = 0.0
slope = 0.001
ustrf = 1.0
ndv = 0


def build_model(idx, ws):
    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = [(11.0, 11, 1.0)]
    ts_times = np.arange(0.0, 12.0, 1.0, dtype=float)
    ts_flows = np.array([1000.0] + [float(q) for q in range(1000, -100, -100)])

    # build MODFLOW 6 files
    name = cases[idx]
    dis_type = dis_types[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
    )
    sim.simulation_data.verify_data = False

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="seconds",
        nper=nper,
        perioddata=tdis_rc,
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(sim, print_option="ALL")

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    if dis_type == "dis":
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
            idomain=1,
        )
    elif dis_type == "disv":
        dis = flopy.mf6.ModflowGwfdisv(
            gwf,
            length_units="meters",
            nlay=nlay,
            ncpl=1,
            nvert=4,
            vertices=vertices,
            cell2d=cell2d,
            top=top,
            botm=botm,
            idomain=1,
        )
    else:
        disukwargs = flopy.utils.gridutil.get_disu_kwargs(
            nlay, nrow, ncol, [delr], [delc], top, [botm]
        )
        dis = flopy.mf6.ModflowGwfdisu(
            gwf,
            vertices=vertices,
            cell2d=cell2d,
            **disukwargs,
        )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf)

    # chd files
    # chd data
    if dis_type == "dis":
        spd = [
            [(0, 0, 0), 0.0],
        ]
    elif dis_type == "disv":
        spd = [
            [(0, 0), 0.0],
        ]
    else:
        spd = [
            [(0,), 0.0],
        ]

    chd = flopy.mf6.modflow.ModflowGwfchd(
        gwf, maxbound=1, stress_period_data=spd, pname="chd-1"
    )

    # sfr file
    if dis_type == "dis":
        if "fail" in name:
            cellid = (2, 2, 2)
        elif "none" in name:
            cellid = "none"
        else:
            cellid = (-1, -1, -1)
    elif dis_type == "disv":
        if "fail" in name:
            cellid = (2, 2)
        elif "none" in name:
            cellid = "none"
        else:
            cellid = (-1, -1)
    else:
        if "fail" in name:
            cellid = (2,)
        elif "none" in name:
            cellid = "none"
        else:
            cellid = (-1,)
    packagedata = []
    for irch in range(nreaches):
        nconn = 1
        if 0 < irch < nreaches - 1:
            nconn += 1
        rp = [
            irch,
            cellid,
            rlen,
            rwid,
            slope,
            top,
            rbth,
            rhk,
            roughness,
            nconn,
            ustrf,
            ndv,
        ]
        packagedata.append(rp)

    if not str(ws).endswith("mf6"):
        packagedata = packagedata[::-1]

    connectiondata = []
    inflow_loc = 0
    ioutflow_loc = nreaches - 1
    for irch in range(nreaches):
        rc = [irch]
        if irch > 0:
            rc.append(irch - 1)
        if irch < nreaches - 1:
            rc.append(-(irch + 1))
        connectiondata.append(rc)

    ts_names = ["inflow"]
    perioddata = [
        [inflow_loc, "inflow", "inflow"],
    ]
    ts_methods = ["linearend"] * len(ts_names)
    ts_data = []
    for t, q in zip(ts_times, ts_flows):
        ts_data.append((t, q))

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        mover=True,
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname="sfr-1",
    )
    fname = f"{name}.sfr.ts"
    sfr.ts.initialize(
        filename=fname,
        timeseries=ts_data,
        time_series_namerecord=ts_names,
        interpolation_methodrecord=ts_methods,
    )
    fname = f"{name}.sfr.obs"
    sfr_obs = {
        f"{fname}.csv": [
            ("inflow", "ext-inflow", (inflow_loc,)),
            ("outflow", "ext-outflow", (ioutflow_loc,)),
        ]
    }
    sfr.obs.initialize(filename=fname, print_input=True, continuous=sfr_obs)

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
    return sim, None


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        compare=None,
        xfail="fail" in name,
    )
    test.run()

import os

import flopy
import numpy as np
import pytest
from cross_section_functions import get_depths
from framework import TestFramework

paktest = "sfr"
cases = [
    "sfr_npt01a",
    "sfr_npt01b",
    "sfr_npt01c",
    "sfr_npt01d",
    "sfr_npt01e",
    "sfr_npt01f",
    "sfr_npt01g",
    "sfr_npt01h",
    "sfr_npt01i",
    "sfr_npt01j",
    "sfr_npt01k",
]

xsect_types = (
    "wide",
    "wide_npoint",
    "slope",
    "slope_r",
    "rectangular",
    "trapezoidal",
    "v",
    "w",
    "v_invalid",
    "|/",
    r"\|",
)

# spatial discretization data
nlay, nrow, ncol = 1, 1, 1
delr, delc = 100.0, 100.0
top = 0.0
botm = -10.0
strt = 0.0

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

np_data = {
    xsect_types[0]: {
        "x": np.array([rwid], dtype=float),
        "h": np.array([0.0], dtype=float),
        "n": np.array([roughness], dtype=float),
    },
    xsect_types[1]: {
        "x": np.array([0.0, rwid], dtype=float),
        "h": np.array([0.0, 0.0], dtype=float),
        "n": np.array([roughness] * 2, dtype=float),
    },
    xsect_types[2]: {
        "x": np.array([0.0, rwid], dtype=float),
        "h": np.array([1.0, 0.0], dtype=float),
        "n": np.array([roughness] * 2, dtype=float),
    },
    xsect_types[3]: {
        "x": np.array([0.0, rwid], dtype=float),
        "h": np.array([0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 2, dtype=float),
    },
    xsect_types[4]: {
        "x": np.array([0.0, 0.0, rwid, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 4, dtype=float),
    },
    xsect_types[5]: {
        "x": np.array([0.0, 0.4 * rwid, 0.6 * rwid, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 4, dtype=float),
    },
    xsect_types[6]: {
        "x": np.array([0.0, 0.5 * rwid, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 3, dtype=float),
    },
    xsect_types[7]: {
        "x": np.array([0.0, 0.2 * rwid, 0.5 * rwid, 0.7 * rwid, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 0.5, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 5, dtype=float),
    },
    xsect_types[8]: {
        "x": np.array([0.0, 0.1 * rwid, 0.5 * rwid, 0.9 * rwid, rwid], dtype=float),
        "h": np.array([1.0, 1.0, 0.0, 1.0, 1.0], dtype=float),
        "n": np.array([roughness] * 5, dtype=float),
    },
    xsect_types[9]: {
        "x": np.array([0.0, 0.0, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 3, dtype=float),
    },
    xsect_types[10]: {
        "x": np.array([0.0, rwid, rwid], dtype=float),
        "h": np.array([1.0, 0.0, 1.0], dtype=float),
        "n": np.array([roughness] * 3, dtype=float),
    },
}


def build_models(idx, test):
    xsect_type = xsect_types[idx]

    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = [(11.0, 11, 1.0)]
    ts_times = np.arange(0.0, 12.0, 1.0, dtype=float)
    ts_flows = np.array([1000.0] + [float(q) for q in range(1000, -100, -100)])

    # build MODFLOW 6 files
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=test.workspace,
    )
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
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
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

    # sfr file
    packagedata = []
    for irch in range(nreaches):
        nconn = 1
        if 0 < irch < nreaches - 1:
            nconn += 1
        rp = [
            irch,
            "none",
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

    sfr_tab = f"{name}.sfr.tab"
    if xsect_type != "wide":
        print("adding n-point cross-section data")
        crosssections = []
        for n in range(nreaches):
            crosssections.append([n, sfr_tab])
    else:
        crosssections = None

    connectiondata = []
    for irch in range(nreaches):
        rc = [irch]
        if irch > 0:
            rc.append(irch - 1)
        if irch < nreaches - 1:
            rc.append(-(irch + 1))
        connectiondata.append(rc)

    ts_names = ["inflow"]
    perioddata = [
        [0, "inflow", "inflow"],
    ]
    ts_methods = ["linearend"] * len(ts_names)
    ts_data = []
    for t, q in zip(ts_times, ts_flows):
        ts_data.append((t, q))

    budpth = f"{name}.{paktest}.cbc"
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        budget_filerecord=budpth,
        mover=True,
        nreaches=nreaches,
        packagedata=packagedata,
        crosssections=crosssections,
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
            ("inflow", "ext-inflow", (0,)),
            ("outflow", "ext-outflow", (nreaches - 1,)),
            ("stage", "stage", (nreaches - 1,)),
            ("depth", "depth", (nreaches - 1,)),
            ("wp", "wet-perimeter", (nreaches - 1,)),
            ("area", "wet-area", (nreaches - 1,)),
        ]
    }
    sfr.obs.initialize(filename=fname, digits=25, print_input=True, continuous=sfr_obs)
    if crosssections is not None:
        stations = np_data[xsect_type]["x"] / rwid
        heights = np_data[xsect_type]["h"]
        roughnesses = np_data[xsect_type]["n"] / roughness
        table = [[x, h, r] for x, h, r in zip(stations, heights, roughnesses)]
        flopy.mf6.ModflowUtlsfrtab(
            gwf,
            nrow=stations.shape[0],
            ncol=3,
            table=table,
            filename=sfr_tab,
        )

    # output control
    budpth = f"{name}.cbc"
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=budpth,
        printrecord=[
            ("BUDGET", "ALL"),
        ],
        saverecord=[
            ("BUDGET", "ALL"),
        ],
    )

    return sim, None


def check_output(idx, test):
    obs_pth = os.path.join(test.workspace, f"{test.name}.sfr.obs.csv")
    obs = flopy.utils.Mf6Obs(obs_pth).get_data()

    assert np.allclose(obs["INFLOW"], np.abs(obs["OUTFLOW"])), (
        "inflow not equal to outflow"
    )

    xs_type = xsect_types[idx]
    xs_d = np_data[xs_type]

    d = get_depths(
        obs["INFLOW"], xs_d["x"], xs_d["h"], roughness=xs_d["n"], slope=slope
    )

    assert np.allclose(obs["DEPTH"], d), (
        "sfr depth not equal to calculated depth.\n"
        f"depth: {obs['DEPTH']}\n"
        f"calculated depth: {d}"
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

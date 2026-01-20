import os

import flopy
import numpy as np
import pytest
from flopy.utils.compare import eval_bud_diff
from framework import TestFramework

paktest = "lak"
budtol = 1e-2
cases = ["ts_lak01"]

# static model data
# spatial discretization
nlay, nrow, ncol = 5, 17, 17
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr = delc = [
    250.00,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    500.00,
    500.00,
    500.00,
    500.00,
    500.00,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    1000.0,
    250.00,
]
top = 500.0
botm = [107.0, 97.0, 87.0, 77.0, 67.0]
idomain = np.ones(shape3d, dtype=int)
idomain[0, 6:11, 6:11] = 0
idomain[1, 7:10, 7:10] = 0

# temporal discretization
nper = 10
sim_time = 15000.0
pertime = sim_time / float(nper)
period_data = []
for n in range(nper):
    period_data.append((pertime, 10, 1.0))

strt = 115.0
icelltype = iconvert = 1
kh, kv, sy, ss = 30.0, [1179.0, 30.0, 30.0, 30.0], 3e-4, 0.2
storage_coefficient = True

# chd data
chd_spd = []
chd_arr = np.linspace(160, 140, ncol)
for k in range(nlay):
    for j in range(nrow):
        if j > 0 and j < nrow - 1:
            chd_spd.append([(k, j, 0), chd_arr[0]])
            chd_spd.append([(k, j, ncol - 1), chd_arr[-1]])
        else:
            for i in range(ncol):
                chd_spd.append([(k, j, i), chd_arr[i]])

# recharge data
recharge = 0.116e-1

# lake data
stage, temp, conc = 110.0, 75.0, 0.5
packagedata = [(0, stage, 57, temp, conc)]
outlets = [(0, 0, -1, "SPECIFIED", -999, -999, -999, -999)]
nlakes = len(packagedata)
noutlets = len(outlets)
connectiondata = [
    (0, 0, (0, 6, 5), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 1, (0, 7, 5), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 2, (0, 8, 5), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 3, (0, 9, 5), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 4, (0, 10, 5), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 5, (0, 5, 6), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 6, (1, 6, 6), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 7, (1, 7, 6), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 8, (1, 7, 6), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 9, (1, 8, 6), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 10, (1, 8, 6), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 11, (1, 9, 6), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 12, (1, 9, 6), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 13, (1, 10, 6), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 14, (0, 11, 6), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 15, (0, 5, 7), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 16, (1, 6, 7), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 17, (1, 6, 7), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 18, (2, 7, 7), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 19, (2, 8, 7), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 20, (2, 9, 7), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 21, (1, 10, 7), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 22, (1, 10, 7), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 23, (0, 11, 7), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 24, (0, 5, 8), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 25, (1, 6, 8), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 26, (1, 6, 8), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 27, (2, 7, 8), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 28, (2, 8, 8), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 29, (2, 9, 8), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 30, (1, 10, 8), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 31, (1, 10, 8), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 32, (0, 11, 8), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 33, (0, 5, 9), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 34, (1, 6, 9), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 35, (1, 6, 9), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 36, (2, 7, 9), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 37, (2, 8, 9), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 38, (2, 9, 9), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 39, (1, 10, 9), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 40, (1, 10, 9), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 41, (0, 11, 9), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 42, (0, 5, 10), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 43, (1, 6, 10), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 44, (1, 7, 10), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 45, (1, 7, 10), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 46, (1, 8, 10), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 47, (1, 8, 10), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 48, (1, 9, 10), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 49, (1, 9, 10), "HORIZONTAL", 0.1, 0, 0, 250, 500),
    (0, 50, (1, 10, 10), "VERTICAL", 0.1, 0, 0, 0, 0),
    (0, 51, (0, 11, 10), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 52, (0, 6, 11), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 53, (0, 7, 11), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 54, (0, 8, 11), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 55, (0, 9, 11), "HORIZONTAL", 0.1, 0, 0, 500, 500),
    (0, 56, (0, 10, 11), "HORIZONTAL", 0.1, 0, 0, 500, 500),
]

stage, evap, runoff, withdrawal, rate = (110.0, 0.0103, 1000.0, 10000.0, -225000.0)
lakeperioddata0 = [
    (0, "status", "active"),
    (0, "stage", stage),
    (0, "rainfall", recharge),
    (0, "evaporation", evap),
    (0, "runoff", runoff),
    (0, "withdrawal", withdrawal),
    (0, "rate", rate),
    (0, "AUXILIARY", "temperature", temp),
    (0, "AUXILIARY", "salinity", conc),
]
lakeperioddata1 = [(0, "rainfall", 0.0)]
lakeperioddata2 = [(0, "rainfall", recharge)]
lakeperioddata3 = [(0, "withdrawal", -rate), (0, "rate", -withdrawal)]
lakeperioddata = {
    0: lakeperioddata0,
    1: lakeperioddata1,
    2: lakeperioddata2,
    3: lakeperioddata3,
}

lakeperioddatats0 = [
    (0, "status", "active"),
    (0, "stage", "stage"),
    (0, "rainfall", "rainfall"),
    (0, "evaporation", "evap"),
    (0, "runoff", "runoff"),
    (0, "withdrawal", "withdrawal"),
    (0, "rate", "outlet"),
    (0, "AUXILIARY", "salinity", "concentration"),
    (0, "AUXILIARY", "temperature", "temperature"),
]
lakeperioddatats1 = [(0, "rainfall", 0.0), (0, "rate", rate)]
lakeperioddatats2 = [(0, "rainfall", "rainfall2"), (0, "rate", "outlet")]
lakeperioddatats3 = [
    (0, "stage", "stage"),
    (0, "rainfall", "rainfall"),
    (0, "evaporation", "evap"),
    (0, "runoff", "runoff"),
    (0, "withdrawal", "outlet2"),
    (0, "rate", "withdrawal2"),
]
lakeperioddatats = {
    0: lakeperioddatats0,
    1: lakeperioddatats1,
    2: lakeperioddatats2,
    3: lakeperioddatats3,
}

ts_names = [
    "stage",
    "rainfall",
    "evap",
    "runoff",
    "withdrawal",
    "outlet",
    "concentration",
    "temperature",
    "rainfall2",
    "outlet2",
    "withdrawal2",
]
ts_methods = ["linearend"] * len(ts_names)

ts_data = []
ts_times = np.arange(0.0, sim_time + 2.0 * pertime, pertime, dtype=float)
for t in ts_times:
    ts_data.append(
        (
            t,
            stage,
            recharge,
            evap,
            runoff,
            withdrawal,
            rate,
            temp,
            conc,
            recharge,
            -rate,
            -withdrawal,
        )
    )

lak_obs = {"lak_obs.csv": [("lake1", "STAGE", (0,))]}


def get_model(ws, name, timeseries=False):
    hdsfile = f"{name}.hds"

    # build the model
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=period_data)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="NONE",
        linear_acceleration="CG",
        outer_maximum=500,
        inner_maximum=100,
        outer_dvclose=1e-6,
        inner_dvclose=1e-3,
        rcloserecord=[0.01, "strict"],
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=kh, icelltype=icelltype)
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        storagecoefficient=storage_coefficient,
        sy=sy,
        ss=ss,
        transient={0: True},
        iconvert=iconvert,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    if timeseries:
        lakpd = lakeperioddatats
    else:
        lakpd = lakeperioddata
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=nlakes,
        noutlets=noutlets,
        print_input=True,
        print_stage=True,
        print_flows=True,
        auxiliary=["temperature", "salinity"],
        packagedata=packagedata,
        connectiondata=connectiondata,
        outlets=outlets,
        perioddata=lakpd,
        pname="lak-1",
    )
    lak.obs.initialize(
        filename=f"{name}.lak.obs",
        digits=20,
        print_input=True,
        continuous=lak_obs,
    )
    if timeseries:
        fname = f"{name}.lak.ts"
        lak.ts.initialize(
            filename=fname,
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=ts_methods,
        )

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "LAST")],
    )

    return sim


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = get_model(ws, name)

    # build MODFLOW 6 files with UZF package
    ws = os.path.join(test.workspace, "mf6")
    mc = get_model(ws, name, timeseries=True)

    return sim, mc


def check_output(idx, test):
    # get ia/ja from binary grid file
    fname = f"{os.path.basename(test.name)}.dis.grb"
    fpth = os.path.join(test.workspace, fname)
    grbobj = flopy.mf6.utils.MfGrdFile(fpth)
    ia = grbobj._datadict["IA"] - 1

    fname = f"{os.path.basename(test.name)}.cbc"

    # open first cbc file
    fpth = os.path.join(test.workspace, fname)
    cobj0 = flopy.utils.CellBudgetFile(fpth, precision="double")

    # open second cbc file
    fpth = os.path.join(test.workspace, "mf6", fname)
    cobj1 = flopy.utils.CellBudgetFile(fpth, precision="double")

    # define file path and evaluate difference
    fname = f"{os.path.basename(test.name)}.cbc.cmp.out"
    fpth = os.path.join(test.workspace, fname)
    eval_bud_diff(fpth, cobj0, cobj1, ia, dtol=0.1)


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

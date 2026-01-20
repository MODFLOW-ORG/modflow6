import os

import flopy
import pytest
from framework import TestFramework

cases = [
    "csub_sub02a",
    "csub_sub02b",
    "csub_sub02c",
    "csub_sub02d",
    "csub_sub02e",
]
cmppth = "mf6_regression"
cg_ske = 1.14e-3 / (500.0 - 20.0)
cg_S = cg_ske * (500.0 - 20.0)
ss = [cg_S, cg_S, cg_ske, cg_ske, cg_S]
storagecoeff = [True, True, False, False, True]
cdelay = [False, True, False, True, True]
ndelaycells = [None, 19, None, 19, 19]

# static model data
nlay, nrow, ncol = 1, 1, 1
nper = 10
perlen = [182.625 for _ in range(nper)]
nstp = [10 for _ in range(nper)]
tsmult = [1.05 for _ in range(nper)]
steady = [False for _ in range(nper)]
delr, delc = 1000.0, 1000.0
top = -100.0
botm = [-600.0]
strt = 0.0
hnoflo = 1e30
hdry = -1e30
hk = 1e6
laytyp = [0]
sy = 0.0

nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-6, 1e-6, 0.97

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

ib = 1

wd = {}
wd6 = {}
for i in range(nper):
    if i % 2 == 0:
        q = -118.3
    else:
        q = 23.66
    d = [[0, 0, 0, q]]
    d6 = [[(0, 0, 0), q]]
    wd[i] = d
    wd6[i] = d6

# sub data
cc = 0.005
cr = 5e-5
void = 0.82
theta = void / (1.0 + void)
kv = 9.72e-6
sgm = 0.0
sgs = 0.0
ini_stress = 0.0
thick = [20.0]
sfe = cr * thick[0]
sfv = cc * thick[0]
lnd = [0]
ldnd = [0]
dp = [[kv, cr, cc]]


def get_model(idx, ws):
    name = cases[idx]
    ss = 1.14e-3
    sc6 = True
    if not storagecoeff[idx]:
        ss /= top - botm[0]
        sc6 = None

    if cdelay[idx]:
        cdelays = "delay"
    else:
        cdelays = "nodelay"

    sub6 = [
        [
            0,
            (0, 0, 0),
            cdelays,
            ini_stress,
            thick[0],
            1.0,
            cc,
            cr,
            theta,
            kv,
            0.0,
            "db01",
        ]
    ]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=f"{name}.nam")

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        storagecoefficient=sc6,
        transient={0: True},
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=1,
        stress_period_data=wd6,
        save_flows=False,
    )

    # csub files
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        print_input=True,
        boundnames=True,
        head_based=True,
        ndelaycells=ndelaycells[idx],
        ninterbeds=1,
        beta=0.0,
        cg_ske_cr=cg_ske,
        packagedata=sub6,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    ws = os.path.join(test.workspace, cmppth)
    mc = get_model(idx, ws)
    return sim, mc


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
        compare="mf6_regression",
    )
    test.run()

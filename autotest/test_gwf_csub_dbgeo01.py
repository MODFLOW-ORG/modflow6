import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_dbgeo01a"]
ndcell = [19]
strt = [0.0]
chdh = [0]
gso = [True]
bso = [True]
# sstate = [None, True]

# comparison data
compdataa = [
    0.6829965295,
    1.3552260620,
    2.0163615680,
    2.6662194730,
    3.3047615550,
    3.9320942150,
    4.5484650870,
    5.1542570520,
    5.7499797990,
    6.3362591640,
    6.9138245770,
    7.4834950010,
    8.0461637950,
    8.6027829980,
    9.1543475090,
    9.7018796520,
    10.2464145700,
    10.7889868300,
    11.3306186100,
    11.8723096200,
    12.4150289800,
    12.9597090600,
    13.5072412800,
    14.0584737000,
    14.6142102900,
    15.1752116300,
    15.7421966800,
    16.3158455400,
    16.8968027900,
    17.4856812200,
    18.0830657300,
    18.6895172900,
    19.3055766500,
    19.9317679500,
    20.5686019500,
    21.2165789800,
    21.8761915600,
    22.5479266700,
    23.2322678200,
    23.9296967200,
    24.6406948800,
    25.3657448400,
    26.1053314000,
    26.8599425100,
    27.6300702200,
    28.4162113200,
    29.2188679600,
    30.0385480600,
    30.8757655600,
    31.7310404400,
    32.6048985000,
    33.4978707400,
    34.4104923500,
    35.3433011400,
    36.2968353000,
    37.2716304000,
    38.2682154600,
    39.2871078800,
    40.3288072900,
    41.3937880000,
    42.4824901100,
    43.5953091600,
    44.7325843800,
    45.8945855500,
    47.0814986400,
    48.2934103700,
    49.5302919900,
    50.7919824900,
    52.0781717600,
    53.3883838800,
    54.7219611100,
    56.0780489100,
    57.4555824900,
    58.8532751700,
    60.2696090000,
    61.7028279600,
    63.1509339500,
    64.6116859300,
    66.0826023000,
    67.5609666300,
    69.0438370300,
    70.5280589700,
    72.0102819100,
    73.4869793800,
    74.9544728300,
    76.4089589400,
    77.8465405400,
    79.2632606700,
    80.6551399000,
    82.0182164100,
    83.3485884600,
    84.6424588900,
    85.8961809100,
    87.1063045100,
    88.2696228000,
    89.3832170600,
    90.4444998100,
    91.4512546100,
    92.4016715000,
    93.2943770700,
]

compdata = {0: compdataa}

# static model data
# spatial discretization
nlay, nrow, ncol = 1, 1, 3
delr, delc = 1.0, 1.0
top = 0.0
bots = [-100.0]
botm = [top] + bots

# temporal discretization
nper = 1
perlen = [1000.0 for _ in range(nper)]
nstp = [100 for _ in range(nper)]
tsmult = [1.05 for _ in range(nper)]
steady = [False for _ in range(nper)]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

hnoflo = 1e30
hdry = -1e30

# idomain data
ib = 1

# npf and sto data
hk = 1e6
laytyp = [1]
ss = 0.0
sy = 0.2

# solver data
nouter, ninner = 50, 100
hclose, rclose, relax = 1e-6, 1e-6, 0.97

# sub data
cc = 100.0
cr = 1.0
void = 0.82
theta = void / (1.0 + void)
kv = 0.025
sgm = 1.7
sgs = 2.2
ini_stress = 1.0
thick = 1.0


# calculate geostatic and effective stress
def calc_stress(sgm0, sgs0, h, bt):
    geo = []
    for k in range(nlay):
        top = bt[k]
        bot = bt[k + 1]
        ht = h
        if ht > top:
            gs = (top - bot) * sgs0
        elif ht < bot:
            gs = (top - bot) * sgm0
        else:
            gs = ((top - ht) * sgm0) + ((ht - bot) * sgs0)
        geo.append(gs)
    # calculate total geostatic stress at bottom of layer
    for k in range(1, nlay):
        geo[k] += geo[k - 1]
    # calculate effective stress at the bottom of the layer
    es = []
    for k in range(nlay):
        es.append(geo[k] - (h - bt[k + 1]))
    return geo, es


def build_models(idx, test):
    c6 = []
    for j in range(0, ncol, 2):
        c6.append([(0, 0, j), chdh[idx]])
    cd6 = {0: c6}

    geo, es = calc_stress(sgm, sgs, strt[idx], botm)
    sub6 = [[0, (0, 0, 1), "delay", -1.0, thick, 1.0, cc, cr, theta, kv, 1.0]]

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    # create iterative model solution and register the gwf model with it
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
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=bots,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt[idx], filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=ss,
        sy=sy,
        storagecoefficient=True,
        transient={0: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=len(c6), stress_period_data=cd6, save_flows=False
    )

    # csub files
    opth = f"{name}.csub.obs"
    ibcsv = f"{name}.ib.strain.csv"
    skcsv = f"{name}.sk.strain.csv"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        print_input=True,
        ndelaycells=ndcell[idx],
        strainib_filerecord=ibcsv,
        straincg_filerecord=skcsv,
        effective_stress_lag=True,
        # compression_indices=True,
        ninterbeds=1,
        sgs=sgs,
        sgm=sgm,
        packagedata=sub6,
        beta=0.0,
        cg_ske_cr=0.0,
    )
    orecarray = {}
    orecarray["csub_obs.csv"] = [
        ("tcomp", "interbed-compaction", (0, 0, 1)),
        ("gs", "gstress-cell", (0, 0, 1)),
        ("es", "estress-cell", (0, 0, 1)),
        ("pcs", "delay-preconstress", (0, 0)),
        ("sk", "sk", (0, 0, 1)),
    ]
    csub_obs_package = csub.obs.initialize(
        filename=opth, digits=10, print_input=True, continuous=orecarray
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    mc = None

    return sim, mc


def check_output(idx, test):
    # MODFLOW 6 total compaction results
    fpth = os.path.join(test.workspace, "csub_obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # set comparison data
    tc0 = compdata[0]

    # calculate maximum absolute error
    diff = tc["TCOMP"] - tc0[:]
    diffmax = np.abs(diff).max()
    dtol = 1e-6
    msg = f"maximum absolute total-compaction difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.comp.cmp.out")
    with open(fpth, "w") as f:
        line = f"{'TOTIM':>15s}"
        line += f" {'CSUB':>15s}"
        line += f" {'MF':>15s}"
        line += f" {'DIFF':>15s}"
        f.write(line + "\n")
        for i in range(diff.shape[0]):
            line = f"{tc0[i]:15g}"
            line += f" {tc['TCOMP'][i]:15g}"
            line += f" {tc0[i]:15g}"
            line += f" {diff[i]:15g}"
            f.write(line + "\n")

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)


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

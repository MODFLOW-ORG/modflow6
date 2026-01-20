import os

import flopy
import numpy as np
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = ["csub_subwt02a", "csub_subwt02b", "csub_subwt02c", "csub_subwt02d"]
timeseries = [True, False, True, False]
cmppth = "mf6_regression"
htol = [None, None, None, None]
dtol = 1e-3
budtol = 1e-2
paktest = "csub"
ump = [None, True, None, True]
ivoid = [0, 1, 0, 1]
gs0 = [0.0, 0.0, 1700.0, 1700.0]

# static model data
pth = str(project_root_path / "autotest" / "data" / "ibc01_ibound.ref")
ib0 = np.genfromtxt(pth)

# temporal discretization
nper = 3
perlen = [1.0, 21915.0, 21915.0]
nstp = [1, 60, 60]
tsmult = [1.0, 1.0, 1.0]
steady = [True, False, False]
ts_times = np.arange(0.0, 60000, 10000.0, dtype=float)

# spatial discretization
nlay, nrow, ncol = 4, ib0.shape[0], ib0.shape[1]
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
nactive = np.count_nonzero(ib0) * nlay

delr, delc = 2000.0, 2000.0
top = 150.0
botm = [50.0, -100.0, -150.0, -350.0]
strt = 100.0

hnoflo = 1e30
hdry = -1e30

# upw data
laytyp = [1, 0, 0, 0]
hk = [4.0, 4.0, 1e-2, 4.0]
k33 = [0.4, 0.4, 1e-2, 0.4]
sy = [0.3, 0.0, 0.0, 0.0]

w1 = [
    (0, 0, 7, 2.2000000e03),
    (0, 1, 4, 2.2000000e03),
    (0, 1, 7, 2.2000000e03),
    (0, 1, 11, 2.2000000e03),
    (0, 2, 3, 2.2000000e03),
    (0, 3, 11, 2.2000000e03),
    (0, 4, 2, 2.2000000e03),
    (0, 4, 12, 2.2000000e03),
    (0, 5, 13, 2.2000000e03),
    (0, 6, 1, 2.2000000e03),
    (0, 13, 1, 2.2000000e03),
    (0, 13, 13, 2.2000000e03),
    (0, 15, 2, 2.2000000e03),
    (0, 15, 12, 2.2000000e03),
    (0, 16, 12, 2.2000000e03),
    (0, 17, 3, 2.2000000e03),
    (0, 17, 11, 2.2000000e03),
    (0, 18, 6, 2.2000000e03),
]
w2 = [
    (0, 0, 7, 2.2000000e03),
    (0, 1, 4, 2.2000000e03),
    (0, 1, 7, 2.2000000e03),
    (0, 1, 11, 2.2000000e03),
    (0, 2, 3, 2.2000000e03),
    (0, 3, 11, 2.2000000e03),
    (0, 4, 2, 2.2000000e03),
    (0, 4, 12, 2.2000000e03),
    (0, 5, 13, 2.2000000e03),
    (0, 6, 1, 2.2000000e03),
    (0, 13, 1, 2.2000000e03),
    (0, 13, 13, 2.2000000e03),
    (0, 15, 2, 2.2000000e03),
    (0, 15, 12, 2.2000000e03),
    (0, 16, 12, 2.2000000e03),
    (0, 17, 3, 2.2000000e03),
    (0, 17, 11, 2.2000000e03),
    (0, 18, 6, 2.2000000e03),
    (1, 8, 9, -7.2000000e04),
    (3, 11, 6, -7.2000000e04),
]
wd = {0: w1, 1: w2, 2: w1}

ws2 = [((1, 8, 9), -7.2000000e04), ((3, 11, 6), -7.2000000e04)]
ws3 = [((1, 8, 9), 0), ((3, 11, 6), 0)]
wd6 = {1: ws2, 2: ws3}

rch0 = [
    ((0, 0, 7), 0.00055),
    ((0, 1, 4), 0.00055),
    ((0, 1, 7), 0.00055),
    ((0, 1, 11), 0.00055),
    ((0, 2, 3), 0.00055),
    ((0, 3, 11), 0.00055),
    ((0, 4, 2), 0.00055),
    ((0, 4, 12), 0.00055),
    ((0, 5, 13), 0.00055),
    ((0, 6, 1), 0.00055),
    ((0, 13, 1), 0.00055),
    ((0, 13, 13), 0.00055),
    ((0, 15, 2), 0.00055),
    ((0, 15, 12), 0.00055),
    ((0, 16, 12), 0.00055),
    ((0, 17, 3), 0.00055),
    ((0, 17, 11), 0.00055),
    ((0, 18, 6), 0.00055),
]
rch6 = {0: rch0}

chd1 = [
    (0, 19, 7, 100.00000, 100.00000),
    (0, 19, 8, 100.00000, 100.00000),
    (1, 19, 7, 100.00000, 100.00000),
    (1, 19, 8, 100.00000, 100.00000),
    (2, 19, 7, 100.00000, 100.00000),
    (2, 19, 8, 100.00000, 100.00000),
    (3, 19, 7, 100.00000, 100.00000),
    (3, 19, 8, 100.00000, 100.00000),
]
cd = {0: chd1}

chd6 = [
    ((0, 19, 7), 100.00000),
    ((0, 19, 8), 100.00000),
    ((1, 19, 7), 100.00000),
    ((1, 19, 8), 100.00000),
    ((2, 19, 7), 100.00000),
    ((2, 19, 8), 100.00000),
    ((3, 19, 7), 100.00000),
    ((3, 19, 8), 100.00000),
]
cd6 = {0: chd6}

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 0.01, 0.97
fluxtol = nactive * rclose

tdis_rc = []
for idx in range(nper):
    tdis_rc.append((perlen[idx], nstp[idx], tsmult[idx]))

# this used to work
# ib = np.zeros((nlay, nrow, ncol), dtype=int)
# for k in range(nlay):
#    ib[k, :, :] = ib0.copy()
ib = []
for k in range(nlay):
    ib.append(ib0.astype(int).copy())

# subwt data
cc = 0.25
cr = 0.01
void = 0.82
theta = void / (1.0 + void)
kv = 999.0
sgm = 1.7
sgs = 2.0
ini_stress = 15.0
delay_flag = 0
thick = [45.0, 70.0, 50.0, 90.0]

zthick = [
    top - botm[0],
    botm[0] - botm[1],
    botm[1] - botm[2],
    botm[2] - botm[3],
]

beta = 0.0
# beta = 4.65120000e-10
gammaw = 9806.65000000
sw = beta * gammaw * theta
ss = [sw for _ in range(nlay)]

swt6 = []
ibcno = 0
for k in range(nlay):
    for i in range(nrow):
        for j in range(ncol):
            iactive = 0
            if ib0[i, j] > 0:
                iactive = 1
            if i == 19 and (j == 7 or j == 8):
                iactive = 0
            if iactive > 0:
                tag = f"{k + 1:02d}_{i + 1:02d}_{j + 1:02d}"
                d = [
                    ibcno,
                    (k, i, j),
                    "nodelay",
                    ini_stress,
                    thick[k],
                    1.0,
                    cc,
                    cr,
                    theta,
                    kv,
                    999.0,
                    tag,
                ]
                swt6.append(d)
                ibcno += 1


def get_model(idx, ws):
    name = cases[idx]

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
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="NEWTON"
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=ib,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=False, icelltype=laytyp, k=hk, k33=k33
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        steady_state={0: True},
        transient={1: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=len(chd6), stress_period_data=cd6, save_flows=False
    )

    # wel files
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(ws2),
        stress_period_data=wd6,
        save_flows=False,
    )
    # recharge file
    flopy.mf6.ModflowGwfrch(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(rch0),
        stress_period_data=rch6,
        save_flows=False,
    )

    # csub files
    gg = []
    if timeseries[idx]:
        sig0v = "geostress"
        ts_methods = ["linearend"]
        ts_data = []
        for t in ts_times:
            ts_data.append((t, gs0[idx]))
    else:
        sig0v = gs0[idx]
    for i in range(nrow):
        for j in range(ncol):
            if ib0[i, j] > 0:
                gg.append([(0, i, j), sig0v])
    sig0 = {0: gg}
    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        # print_input=True,
        # interbed_stress_offset=True,
        boundnames=False,
        compression_indices=True,
        update_material_properties=ump[idx],
        effective_stress_lag=True,
        ninterbeds=len(swt6),
        sgs=sgs,
        sgm=sgm,
        beta=beta,
        gammaw=gammaw,
        cg_ske_cr=0.0,
        cg_theta=theta,
        packagedata=swt6,
        maxsig0=len(gg),
        stress_period_data=sig0,
    )
    if timeseries[idx]:
        fname = f"{name}.csub.ts"
        csub.ts.initialize(
            filename=fname,
            timeseries=ts_data,
            time_series_namerecord=[sig0v],
            interpolation_methodrecord=ts_methods,
        )

    cobs = [
        ("w1l1", "interbed-compaction", (89,)),
        ("w1l2", "interbed-compaction", (299,)),
        ("w1l3", "interbed-compaction", (509,)),
        ("w1l4", "interbed-compaction", (719,)),
        ("w2l1", "interbed-compaction", (130,)),
        ("w2l2", "interbed-compaction", (340,)),
        ("w2l3", "interbed-compaction", (550,)),
        ("w2l4", "interbed-compaction", (760,)),
        ("s1l1", "coarse-compaction", (0, 8, 9)),
        ("s1l2", "coarse-compaction", (1, 8, 9)),
        ("s1l3", "coarse-compaction", (2, 8, 9)),
        ("s1l4", "coarse-compaction", (3, 8, 9)),
        ("s2l1", "coarse-compaction", (0, 11, 6)),
        ("s2l2", "coarse-compaction", (1, 11, 6)),
        ("s2l3", "coarse-compaction", (2, 11, 6)),
        ("s2l4", "coarse-compaction", (3, 11, 6)),
        ("c1l1", "compaction-cell", (0, 8, 9)),
        ("c1l2", "compaction-cell", (1, 8, 9)),
        ("c1l3", "compaction-cell", (2, 8, 9)),
        ("c1l4", "compaction-cell", (3, 8, 9)),
        ("c2l1", "compaction-cell", (0, 11, 6)),
        ("c2l2", "compaction-cell", (1, 11, 6)),
        ("c2l3", "compaction-cell", (2, 11, 6)),
        ("c2l4", "compaction-cell", (3, 11, 6)),
        ("w2l4q", "csub-cell", (3, 11, 6)),
        ("gs1", "gstress-cell", (0, 8, 9)),
        ("es1", "estress-cell", (0, 8, 9)),
        ("pc1", "preconstress-cell", (0, 8, 9)),
        ("gs2", "gstress-cell", (1, 8, 9)),
        ("es2", "estress-cell", (1, 8, 9)),
        ("pc2", "preconstress-cell", (1, 8, 9)),
        ("gs3", "gstress-cell", (2, 8, 9)),
        ("es3", "estress-cell", (2, 8, 9)),
        ("pc3", "preconstress-cell", (2, 8, 9)),
        ("gs4", "gstress-cell", (3, 8, 9)),
        ("es4", "estress-cell", (3, 8, 9)),
        ("pc4", "preconstress-cell", (3, 8, 9)),
        ("sk1l2", "ske-cell", (1, 8, 9)),
        ("sk2l4", "ske-cell", (3, 11, 6)),
        ("t1l2", "theta", (1, 8, 9)),
    ]

    orecarray = {"csub_obs.csv": cobs}
    csub_obs_package = csub.obs.initialize(
        filename=opth, digits=10, print_input=True, continuous=orecarray
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    sim = get_model(idx, ws)

    # build comparison files
    ws = os.path.join(test.workspace, cmppth)
    mc = get_model(idx, ws)

    return sim, mc


def check_output(idx, test):
    # MODFLOW 6 total compaction results
    fpth = os.path.join(test.workspace, "csub_obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # comparison total compaction results
    cpth = cmppth
    fpth = os.path.join(test.workspace, cmppth, "csub_obs.csv")
    try:
        tc0 = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    loctag = "W2L4"
    diff = tc[loctag] - tc0[loctag]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-compaction difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.comp.cmp.out")
    f = open(fpth, "w")
    line = f"{'TOTIM':>15s}"
    line += f" {'CSUB':>15s}"
    line += f" {'MF':>15s}"
    line += f" {'DIFF':>15s}"
    f.write(line + "\n")
    for i in range(diff.shape[0]):
        line = f"{tc0['time'][i]:15g}"
        line += f" {tc[loctag][i]:15g}"
        line += f" {tc0[loctag][i]:15g}"
        line += f" {diff[i]:15g}"
        f.write(line + "\n")
    f.close()

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)

    # compare budgets
    cbc_compare(test)

    return


# compare cbc and lst budgets
def cbc_compare(test):
    # open cbc file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.cbc")
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    # build list of cbc data to retrieve
    avail = cobj.get_unique_record_names()
    cbc_bud = []
    bud_lst = []
    for t in avail:
        if isinstance(t, bytes):
            t = t.decode()
        t = t.strip()
        if paktest in t.lower():
            cbc_bud.append(t)
            bud_lst.append(f"{t}_IN")
            bud_lst.append(f"{t}_OUT")

    # get results from listing file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.lst")
    budl = flopy.utils.Mf6ListBudget(fpth)
    names = list(bud_lst)
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]
    d = np.recarray(nbud, dtype=dtype)
    for key in bud_lst:
        d[key] = 0.0

    # get data from cbc dile
    kk = cobj.get_kstpkper()
    times = cobj.get_times()
    for i, (k, t) in enumerate(zip(kk, times)):
        for text in cbc_bud:
            qin = 0.0
            qout = 0.0
            v = cobj.get_data(kstpkper=k, text=text)[0]
            if isinstance(v, np.recarray):
                vt = np.zeros(size3d, dtype=float)
                for jdx, node in enumerate(v["node"]):
                    vt[node - 1] += v["q"][jdx]
                v = vt.reshape(shape3d)
            for kk in range(v.shape[0]):
                for ii in range(v.shape[1]):
                    for jj in range(v.shape[2]):
                        vv = v[kk, ii, jj]
                        if vv < 0.0:
                            qout -= vv
                        else:
                            qin += vv
            d["totim"][i] = t
            d["time_step"][i] = k[0]
            d["stress_period"] = k[1]
            key = f"{text}_IN"
            d[key][i] = qin
            key = f"{text}_OUT"
            d[key][i] = qout

    diff = np.zeros((nbud, len(bud_lst)), dtype=float)
    for i, key in enumerate(bud_lst):
        diff[:, i] = d0[key] - d[key]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-budget difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.bud.cmp.out")
    with open(fpth, "w") as f:
        for i in range(diff.shape[0]):
            if i == 0:
                line = f"{'TIME':>10s}"
                for key in bud_lst:
                    line += f"{key + '_LST':>25s}"
                    line += f"{key + '_CBC':>25s}"
                    line += f"{key + '_DIF':>25s}"
                f.write(line + "\n")
            line = f"{d['totim'][i]:10g}"
            for ii, key in enumerate(bud_lst):
                line += f"{d0[key][i]:25g}"
                line += f"{d[key][i]:25g}"
                line += f"{diff[i, ii]:25g}"
            f.write(line + "\n")

    if diffmax > budtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        htol=htol[idx],
        compare="mf6_regression",
    )
    test.run()

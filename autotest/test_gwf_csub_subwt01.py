import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_subwt01a", "csub_subwt01b", "csub_subwt01c", "csub_subwt01d"]
cmppth = "mf6_regression"
htol = [None for _ in cases]
dtol = 1e-3
budtol = 1e-2
paktest = "csub"
ump = [None, True, None, True]
ivoid = [0, 1, 0, 1]
gs0 = [0.0, 0.0, 1700.0, 1700.0]

# temporal discretization
nper = 3
perlen = [1.0, 21915.0, 21915.0]
nstp = [1, 60, 60]
tsmult = [1.0, 1.0, 1.0]
steady = [True, False, False]

# spatial discretization
nlay, nrow, ncol = 1, 1, 3
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol

delr, delc = 2000.0, 2000.0
top = 150.0
botm = [-350.0]
strt = 100.0

hnoflo = 1e30
hdry = -1e30

# upw data
laytyp = [1]
hk = [4.0]
sy = [0.3]

hstart = [100.0, 50.0]

chd0 = [
    (0, 0, 0, hstart[0] + 1, hstart[0] + 1),
    (0, 0, 2, hstart[0], hstart[0]),
]
chd1 = [
    (0, 0, 0, hstart[1] + 1, hstart[1] + 1),
    (0, 0, 2, hstart[1], hstart[1]),
]
cd = {0: chd0, 1: chd1, 2: chd0}

chd6_0 = [((0, 0, 0), hstart[0] + 1), ((0, 0, 2), hstart[0])]
chd6_1 = [((0, 0, 0), hstart[1] + 1), ((0, 0, 2), hstart[1])]
cd6 = {0: chd6_0, 1: chd6_1, 2: chd6_0}

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 0.01, 0.97
fluxtol = rclose

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# this used to work
# ib = np.zeros((nlay, nrow, ncol), dtype=int)
# for k in range(nlay):
ib = [1]

# subwt data
cc = 0.25
cr = 0.01
void = 0.82
theta = void / (1.0 + void)
kv = 999.0
sgm = 1.7
sgs = 2.0
ini_stress = 0.0  # 15.0
delay_flag = 0
thick = [45.0, 70.0, 50.0, 90.0]

zthick = [top - botm[0]]

beta = 0.0
# beta = 4.65120000e-10
gammaw = 9806.65000000
sw = beta * gammaw * theta
ss = [sw for k in range(nlay)]

swt6 = []
ibcno = 0
for k in range(len(thick)):
    for i in range(nrow):
        for j in range(ncol):
            iactive = 0
            if j == 1:
                iactive = 1
            if iactive > 0:
                tag = f"{1:02d}_{i + 1:02d}_{j + 1:02d}"
                d = [
                    ibcno,
                    (0, i, j),
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
        sim_name=name,
        memory_print_option="all",
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True, print_input=True)

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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)
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
        gwf, maxbound=len(chd6_0), stress_period_data=cd6, save_flows=False
    )

    # csub files
    gg = []
    for i in range(nrow):
        for j in range(ncol):
            gg.append([(0, i, j), gs0[idx]])
    sig0 = {0: gg}
    opth = f"{name}.csub.obs"
    fcgstrain = f"{name}.csub.strain.cg.csv"
    fibstrain = f"{name}.csub.strain.ib.csv"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        straincg_filerecord=fcgstrain,
        strainib_filerecord=fibstrain,
        boundnames=True,
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
    orecarray = {}
    orecarray["csub_obs.csv"] = [
        ("w1l1", "compaction-cell", (0, 0, 1)),
        ("w1l1t", "csub-cell", (0, 0, 1)),
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
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    sim = get_model(idx, test.workspace)

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

    # Comparison total compaction results
    cpth = cmppth
    fpth = os.path.join(test.workspace, cpth, "csub_obs.csv")
    try:
        tc0 = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    loctag = "W1L1"
    diff = tc[loctag] - tc0[loctag]
    diffmax = np.abs(diff).max()
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
            line = f"{tc0['time'][i]:15g}"
            line += f" {tc[loctag][i]:15g}"
            line += f" {tc0[loctag][i]:15g}"
            line += f" {diff[i]:15g}"
            f.write(line + "\n")

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)

    # compare budgets
    cbc_compare(test)


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

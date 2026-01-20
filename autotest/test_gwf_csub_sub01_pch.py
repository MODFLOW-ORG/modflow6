import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_sub01_pch"]
paktest = "csub"
budtol = 1e-2
compdir = "mf6"
ndcell = [19] * len(cases)

# static model data
# spatial discretization
nlay, nrow, ncol = 1, 1, 3
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1.0, 1.0
top = 0.0
botm = [-100.0]

# temporal discretization
nper = 1
perlen = [1000.0 for _ in range(nper)]
nstp = [100 for _ in range(nper)]
tsmult = [1.05 for _ in range(nper)]
steady = [False for _ in range(nper)]

strt = 0.0
strt6 = 1.0
hnoflo = 1e30
hdry = -1e30
hk = 1e6
laytyp = [0]
S = 1e-4
sy = 0.0

nouter, ninner = 1000, 500
hclose, rclose, relax = 1e-12, 1e-6, 0.97

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

ib = 1

c = []
c6 = []
for j in range(0, ncol, 2):
    c.append([0, 0, j, strt, strt])
    c6.append([(0, 0, j), strt])
cd = {0: c}
cd6 = {0: c6}

# sub data
ndb = 1
nndb = 0
cc = 1.0e-2
cr = 1.0e-5
theta = 0.45
void = theta / (1.0 - theta)
kv = 2.5e-6
sgm = 0.0
sgs = 0.0
ini_stress = 1.0
ini_head = 1.0
thick = [1.0]


def build_models(idx, test):
    sim = get_model(idx, test.workspace, pch=True)

    # build MODFLOW-6 with constant material properties
    pth = os.path.join(test.workspace, compdir)
    mc = get_model(idx, pth)

    return sim, mc


def get_model(idx, dir, pch=None):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    sim.name_file.memory_print_option = "all"

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
        rcloserecord=f"{rclose} strict",
        linear_acceleration="bicgstab",
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
        storagecoefficient=True,
        transient={0: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=len(c6), stress_period_data=cd6, save_flows=False
    )

    # csub files
    if pch is None:
        sub6 = [
            [
                0,
                (0, 0, 1),
                "delay",
                -ini_stress,
                thick[0],
                1.0,
                cc,
                cr,
                theta,
                kv,
                ini_head,
            ]
        ]
    else:
        sub6 = [
            [
                0,
                (0, 0, 1),
                "delay",
                ini_stress,
                thick[0],
                1.0,
                cc,
                cr,
                theta,
                kv,
                ini_head,
            ]
        ]

    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        initial_preconsolidation_head=pch,
        print_input=True,
        save_flows=True,
        ndelaycells=ndcell[idx],
        ninterbeds=1,
        beta=0.0,
        cg_ske_cr=0.0,
        packagedata=sub6,
    )
    obs = [
        ("tcomp", "compaction-cell", (0, 0, 1)),
        ("thick", "thickness", (0,)),
        ("theta", "theta", (0,)),
        ("sk", "sk", (0, 0, 1)),
    ]
    orecarray = {}
    orecarray["csub_obs.csv"] = obs

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
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def check_output(idx, test):
    # MODFLOW 6 compaction results
    fpth = os.path.join(test.workspace, "csub_obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # MODFLOW 6 base compaction results
    fpth = os.path.join(test.workspace, compdir, "csub_obs.csv")
    try:
        tcb = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    diff = tc["TCOMP"] - tcb["TCOMP"]
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
            line = f"{tc['time'][i]:15g}"
            line += f" {tc['TCOMP'][i]:15g}"
            line += f" {tcb['TCOMP'][i]:15g}"
            line += f" {diff[i]:15g}"
            f.write(line + "\n")

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol:15.7g}"
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

    dtol = 1e-6

    if diffmax > budtol:
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

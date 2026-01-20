import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_sk02a", "csub_sk02b", "csub_sk02c", "csub_sk02d"]
constantcv = [True for _ in range(len(cases))]
cmppths = ["mf6_regression" for _ in range(len(cases))]
tops = [150.0 for _ in range(len(cases))]
newtons = [True for _ in range(len(cases))]
ump = [None, None, True, True]
iump = [0, 0, 1, 1]
eslag = [True for _ in range(len(cases))]
icrcc = [0, 1, 0, 1]
htol = [None for _ in range(len(cases))]
dtol = 1e-3
bud_lst = [
    "CSUB-CGELASTIC_IN",
    "CSUB-CGELASTIC_OUT",
    "CSUB-WATERCOMP_IN",
    "CSUB-WATERCOMP_OUT",
]

# static model data
nlay, nrow, ncol = 3, 10, 10
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
delr, delc = 1000.0, 2000.0
top = 150.0
botm = [-100, -150.0, -350.0]
zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
strt = 100.0
hnoflo = 1e30
hdry = -1e30

# calculate hk
hk1fact = 1.0 / zthick[1]
hk1 = np.ones((nrow, ncol), dtype=float) * 0.5 * hk1fact
hk1[0, :] = 1000.0 * hk1fact
hk1[-1, :] = 1000.0 * hk1fact
hk1[:, 0] = 1000.0 * hk1fact
hk1[:, -1] = 1000.0 * hk1fact
hk = [20.0, hk1, 5.0]

# calculate vka
vka = [1e6, 7.5e-5, 1e6]

# set rest of npf variables
laytyp = [1, 0, 0]
laytypu = [4, 0, 0]
sy = [0.1, 0.0, 0.0]

nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# all cells are active
ib = 1

# chd data
c = []
c6 = []
ccol = [3, 4, 5, 6]
for j in ccol:
    c.append([0, nrow - 1, j, strt, strt])
    c6.append([(0, nrow - 1, j), strt])
cd = {0: c}
cd6 = {0: c6}
maxchd = len(cd[0])

# pumping well data
wr = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3]
wc = [0, 1, 8, 9, 0, 9, 0, 9, 0, 0]
wrp = [2, 2, 3, 3]
wcp = [5, 6, 5, 6]
wq = [-14000.0, -8000.0, -5000.0, -3000.0]
d = []
d6 = []
for r, c, q in zip(wrp, wcp, wq):
    d.append([2, r, c, q])
    d6.append([(2, r, c), q])
wd = {1: d}
wd6 = {1: d6}
maxwel = len(wd[1])

# recharge data
q = 3000.0 / (delr * delc)
v = np.zeros((nrow, ncol), dtype=float)
for r, c in zip(wr, wc):
    v[r, c] = q
rech = {0: v}

# static ibc and sub data
sgm = 1.7
sgs = 2.0
void = 0.82
ini_stress = 15.0
theta = void / (1.0 + void)
sw = 4.65120000e-10 * 9806.65000000 * theta

# no delay bed data
lnd = [0, 1, 2]
thicknd0 = [zthick[0], zthick[1], zthick[2]]
cr = [0.01, 0.005, 0.01]
sske = [6e-6, 3e-6, 6e-6]

# subwt output data
ds16 = [0, 0, 0, 2052, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
ds17 = [
    0,
    10000,
    0,
    10000,
    0,
    0,
    1,
    1,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
]


def get_model(idx, ws):
    name = cases[idx]
    newton = newtons[idx]
    newtonoptions = None
    imsla = "CG"
    if newton:
        newtonoptions = "NEWTON"
        imsla = "BICGSTAB"

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
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
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    if icrcc[idx] == 0:
        sc = cr
        compression_indices = True
    else:
        sc = sske
        compression_indices = None
    # water compressibility cannot be compared for cases where the material
    # properties are adjusted since the porosity changes in mf6
    if iump[idx] == 0:
        beta = 4.6512e-10
        wc = sw
    else:
        beta = 0.0
        wc = 0.0

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions=newtonoptions)

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
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        # dev_modflowusg_upstream_weighted_saturation=True,
        icelltype=laytyp,
        k=hk,
        k33=vka,
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        storagecoefficient=True,
        steady_state={0: True},
        transient={1: True},
    )

    # recharge
    rch = flopy.mf6.ModflowGwfrcha(gwf, readasarrays=True, recharge=rech)

    # wel file
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=maxwel,
        stress_period_data=wd6,
        save_flows=False,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
    )
    # ibc files
    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        update_material_properties=ump[idx],
        effective_stress_lag=eslag[idx],
        save_flows=True,
        ninterbeds=0,
        compression_indices=compression_indices,
        sgm=sgm,
        sgs=sgs,
        cg_theta=theta,
        cg_ske_cr=sc,
        beta=beta,
        packagedata=None,
    )
    orecarray = {}
    orecarray["csub_obs.csv"] = [
        ("tcomp1", "compaction-cell", (0, 4, 4)),
        ("tcomp2", "compaction-cell", (1, 4, 4)),
        ("tcomp3", "compaction-cell", (2, 4, 4)),
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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    sim = get_model(idx, test.workspace)

    # build comparison files
    cpth = cmppths[idx]
    ws = os.path.join(test.workspace, cpth)
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
    cpth = cmppths[idx]
    fpth = os.path.join(test.workspace, cpth, "csub_obs.csv")
    try:
        tc0 = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    diff = tc["TCOMP3"] - tc0["TCOMP3"]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-compaction difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.comp.cmp.out")
    with open(fpth, "w") as f:
        for i in range(diff.shape[0]):
            line = f"{tc0['time'][i]:10.2g}"
            line += f"{tc['TCOMP3'][i]:10.2g}"
            line += f"{tc0['TCOMP3'][i]:10.2g}"
            line += f"{diff[i]:10.2g}"
            f.write(line + "\n")

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)

    # get results from listing file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.lst")
    budl = flopy.utils.Mf6ListBudget(fpth)
    names = list(bud_lst)
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]

    # get results from cbc file
    cbc_bud = ["CSUB-CGELASTIC", "CSUB-WATERCOMP"]
    d = np.recarray(nbud, dtype=dtype)
    for key in bud_lst:
        d[key] = 0.0
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.cbc")
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    kk = cobj.get_kstpkper()
    times = cobj.get_times()
    for i, (k, t) in enumerate(zip(kk, times)):
        for text in cbc_bud:
            qin = 0.0
            qout = 0.0
            v = cobj.get_data(kstpkper=k, text=text)[0]
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
    for idx, key in enumerate(bud_lst):
        diff[:, idx] = d0[key] - d[key]
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
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        htol=htol[idx],
        compare="mf6_regression",
    )
    test.run()

import datetime
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_sk03a"]
constantcv = [True for _ in range(len(cases))]
cmppths = ["mf6_regression" for _ in range(len(cases))]
newtons = [True for _ in range(len(cases))]
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
# temporal discretization
nper = 2
sec2day = 86400.0
day2sec = 1.0 / sec2day
nsec = 33 * 60
perlen = np.array([1.0, nsec])
totim = perlen.sum() - perlen[0]
nstp = [1, nsec * 2]
tsmult = [1.0, 1.00]
steady = [True] + [False for _ in range(nper - 1)]

# spatial discretization
ft2m = 1.0 / 3.28081
nlay, nrow, ncol = 3, 21, 20
delr = np.ones(ncol, dtype=float) * 0.5
for i in range(1, ncol):
    delr[i] = min(delr[i - 1] * 1.2, 15.0)
delc = 50.0
top = 0.0
botm = np.array([-40, -70.0, -100.0], dtype=float) * ft2m
zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
strt = -35.0 * ft2m
hnoflo = 1e30
hdry = -1e30

# calculate hk
hk = np.array([5.0, 0.001, 15.0]) * ft2m * day2sec

# calculate vka
vka = hk.copy() * 0.1

# set rest of npf variables
laytyp = [1, 0, 0]
sy = [0.1, 0.05, 0.25]

nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# all cells are active
ib = 1

# chd data
finish = strt - totim * 0.026 / (60.0 * 60.0 * 3.28081)
c = []
c6 = []
ccol = [ncol - 1]
for k in [0, nlay - 1]:
    for i in range(nrow):
        for j in ccol:
            c.append([k, i, j, strt, strt])
            c6.append([(k, i, j), "chd"])
cd = {0: c}
cd6 = {0: c6}
maxchd = len(cd[0])

# static csub and subwt data
ump = True
sgm = 1.7
sgs = 2.0
ini_stress = 15.0
theta = [0.25, 0.5, 0.30]
beta = 4.65120000e-10

# no delay bed data
lnd = [0, 1, 2]
thicknd0 = [zthick[0], zthick[1], zthick[2]]
sske = np.array([1.0e-5, 2.0e-4, 9.0e-8]) / ft2m

gl0 = np.zeros((nrow, ncol), dtype=float)
jj = 0
gl0[0, jj] = 1.86
sig0 = np.zeros((nlay, nrow, ncol), dtype=float)
sig0[0, 0, jj] = 1.86
tsnames = []
sig0 = []
for i in range(nrow):
    tsname = f"TLR{i + 1:02d}"
    tsnames.append(tsname)
    sig0.append([(0, i, 0), tsname])
tsname = "FR"
tsnames.append(tsname)
sig0.append([(0, 9, 0), tsname])

datestart = datetime.datetime.strptime("03/21/1938 00:00:00", "%m/%d/%Y %H:%M:%S")
train1 = 2.9635  # 3.9009
train2 = 2.8274
fcar1 = 0.8165
fcar2 = 1.63293447
fcar3 = 2.8274  # 2.9635 #2.4494
icnt = np.zeros((nrow), dtype=int)
v = []
tstime = []
i0 = 0
dt = 15.0
t = []
ton = 6.0 * 60.0
for i in range(nrow + 1):
    vv = 0.0
    t.append(vv)
tstime.append(ton)
v.append(t)
ton += 1.0
train = train1
fend = train1  # fcar3
fremain = 0.0
while True:
    if i0 < nrow:
        icnt[i0] = 1
    t = []
    for i in range(nrow):
        vv = 0.0
        if icnt[i] > 0 and icnt[i] < 7:
            if i == i0:
                vv = train
            elif i == i0 - 5:
                vv = fend
            else:
                vv = train  # fcar3
        t.append(vv)
    t.append(fremain)
    ton += dt
    tstime.append(ton)
    v.append(t)
    if i0 == 13:
        ton += 525
        tstime.append(ton)
        v.append(t)
    if icnt[9] >= 6:
        train = train2
        fend = train2  # fcar2
        fremain = fcar1
    icnt[icnt > 0] += 1
    i0 += 1
    if icnt[nrow - 1] == 7:
        break
tsv = np.array(v)

# sig0 = [[(0, 0, 0), tsname]]

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
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="SECONDS", nper=nper, perioddata=tdis_rc
    )

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
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    sc = sske
    compression_indices = None

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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=vka,
    )
    # gwf obs
    jj = 10
    obs_recarray = {
        "gwf_obs.csv": [
            ("t1_1_1", "HEAD", (0, 10, 0)),
            ("t2_1_1", "HEAD", (1, 10, 0)),
            ("t3_1_1", "HEAD", (2, 10, 0)),
            ("w1_1_1", "HEAD", (0, 10, jj)),
            ("w2_1_1", "HEAD", (1, 10, jj)),
            ("w3_1_1", "HEAD", (2, 10, jj)),
            ("ICF1_1_1", "FLOW-JA-FACE", (0, 10, 0), (0, 10, 1)),
            ("ICF2_1_1", "FLOW-JA-FACE", (1, 10, 0), (1, 10, 1)),
            ("ICF3_1_1", "FLOW-JA-FACE", (2, 10, 0), (2, 10, 1)),
        ],
        "gwf_calib_obs.csv": [("w3_1_1", "HEAD", (2, 10, jj))],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        gwf,
        pname="head_obs",
        digits=10,
        print_input=True,
        continuous=obs_recarray,
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

    # create chd time series
    chnam = f"{name}.ch.ts"
    chd_ts = [(0.0, strt), (1.0, strt), (perlen.sum(), finish)]

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
    )
    # initialize time series
    chd.ts.initialize(
        filename=chnam,
        timeseries=chd_ts,
        time_series_namerecord=["CHD"],
        interpolation_methodrecord=["linear"],
    )

    # create load time series file with load
    csub_ts = []
    csubnam = f"{name}.load.ts"

    fopth = os.path.join(ws, "ts.csv")
    fo = open(fopth, "w")

    line2 = "TIME"
    for tsname in tsnames:
        line2 += f",{tsname}"
    fo.write(f"{line2}\n")

    d = [0.0]
    line = f" {0.0:12.6e}"
    dateon = datestart + datetime.timedelta(seconds=0)
    datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
    line2 = f"{datestr}"
    for tsname in tsnames:
        d.append(0.0)
        line2 += ",0.0000"
    csub_ts.append(tuple(d))
    fo.write(f"{line2}\n")

    ton = 1.0
    d = [ton]
    dateon = datestart + datetime.timedelta(seconds=ton)
    datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
    line2 = f"{datestr}"
    for tsname in tsnames:
        d.append(0.0)
        line2 += ",0.0000"
    csub_ts.append(tuple(d))
    fo.write(f"{line2}\n")

    for i in range(tsv.shape[0]):
        ton = tstime[i]
        d = [ton]
        dateon = datestart + datetime.timedelta(seconds=ton)
        datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
        line2 = f"{datestr}"
        for j in range(tsv.shape[1]):
            d.append(tsv[i, j])
            line2 += f",{tsv[i, j]:6.4f}"
        csub_ts.append(tuple(d))
        fo.write(f"{line2}\n")

    ton += 1
    d = [ton]
    dateon = datestart + datetime.timedelta(seconds=ton)
    datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
    line2 = f"{datestr}"
    for tsname in tsnames:
        if tsname == "FR":
            d.append(fcar1)
            line2 += f",{fcar1:6.4f}"
        else:
            d.append(0.0)
            line2 += ",0.0000"
    csub_ts.append(tuple(d))
    fo.write(f"{line2}\n")

    ton += (i + 1.0) * 86400.0
    d = [ton]
    dateon = datestart + datetime.timedelta(minutes=29)
    datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
    line2 = f"{datestr}"
    for tsname in tsnames:
        if tsname == "FR":
            d.append(fcar1)
            line2 += f",{fcar1:6.4f}"
        else:
            d.append(0.0)
            line2 += ",0.0000"
    csub_ts.append(tuple(d))
    fo.write(f"{line2}\n")

    ton = 100.0 * sec2day
    d = [ton]
    dateon = datestart + datetime.timedelta(minutes=30)
    datestr = dateon.strftime("%m/%d/%Y %H:%M:%S")
    line2 = f"{datestr}"
    for tsname in tsnames:
        if tsname == "FR":
            d.append(fcar1)
            line2 += f",{fcar1:6.4f}"
        else:
            d.append(fcar1)
            line2 += ",0.0000"
    csub_ts.append(tuple(d))
    fo.write(f"{line2}\n")

    # close ts,csv
    fo.close()

    # csub files
    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        print_input=True,
        update_material_properties=ump,
        effective_stress_lag=True,
        save_flows=True,
        ninterbeds=0,
        maxsig0=len(sig0),
        compression_indices=compression_indices,
        sgm=sgm,
        sgs=sgs,
        cg_theta=theta,
        cg_ske_cr=sc,
        beta=beta,
        packagedata=None,
        stress_period_data={0: sig0},
    )
    # initialize time series
    csub.ts.initialize(
        filename=csubnam,
        timeseries=csub_ts,
        time_series_namerecord=tsnames,
        interpolation_methodrecord=["linear" for t in tsnames],
    )

    # add observations
    orecarray = {}
    jj = 10
    orecarray["csub_obs.csv"] = [
        ("tc1", "compaction-cell", (0, 10, 0)),
        ("tc2", "compaction-cell", (1, 10, 0)),
        ("tc3", "compaction-cell", (2, 10, 0)),
        ("wc1", "compaction-cell", (0, 10, jj)),
        ("wc2", "compaction-cell", (1, 10, jj)),
        ("wc3", "compaction-cell", (2, 10, jj)),
        ("tes3", "estress-cell", (2, 10, 0)),
        ("fes3", "estress-cell", (2, 9, 0)),
        ("tgs3", "gstress-cell", (2, 10, 0)),
        ("fgs3", "gstress-cell", (2, 9, 0)),
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
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


# SUB package problem 3
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

    if diffmax > dtol:
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
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        htol=htol[idx],
        compare="mf6_regression",
    )
    test.run()

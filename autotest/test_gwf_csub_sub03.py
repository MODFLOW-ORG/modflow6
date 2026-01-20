import os

import flopy
import numpy as np
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = ["csub_sub03a", "csub_sub03b"]
cmppth = "mf6_regression"
cvopt = [None, None, None]
constantcv = [True, True]
ndelaybeds = [0, 2]
ndelaycells = [None, 39]
htol = [None, None, None]
dtol = 1e-3
bud_lst = [
    "CSUB-ELASTIC_IN",
    "CSUB-INELASTIC_IN",
    "CSUB-ELASTIC_OUT",
    "CSUB-INELASTIC_OUT",
]

# static model data
nlay, nrow, ncol = 3, 10, 10
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
delr, delc = 1000.0, 2000.0
top = 0.0
botm = [-100, -150.0, -350.0]
zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
strt = 0.0
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
laytyp = [0, 0, 0]  # [1, 0, 0]
ffrac = [0.6, 0.0, 0.6]
sy = [0.1, 0.0, 0.0]
ss = [3e-6, 0.0, 3e-6]

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
sgm = 0.0
sgs = 0.0
omega = 1.0

# no delay bed data
nndb = 3
lnd = [0, 1, 2]
hc = [-7.0, -7.0, -7.0]
thicknd0 = [15.0, 50.0, 30.0]
ccnd0 = [6e-4, 3e-4, 6e-4]
crnd0 = [6e-6, 3e-6, 6e-6]
sfv = []
sfe = []
for k in range(nlay):
    sfv.append(ccnd0[k] * thicknd0[k])
    sfe.append(crnd0[k] * thicknd0[k])

# delay bed data
nmz = 1
ldnd = [0, 2]
kv = 1e-6
void = 0.82
theta = void / (1.0 + void)
cc = 6e-4
cr = 6e-6
dp = [[kv, cr, cc]]
rnb = [7.635, 17.718]
dhc = [-7.0, -7.0]
dz = [5.894, 5.08]
nz = [1, 1]
dstart = []
for k in ldnd:
    pth = str(project_root_path / "autotest" / "data" / f"ibc03_dstart{k + 1}.ref")
    v = np.genfromtxt(pth)
    dstart.append(v.copy())


# SUB package problem 3
def get_model(idx, ws):
    name = cases[idx]
    # ibc packagedata container counter
    sub6 = []
    ibcno = 0

    # create no delay bed packagedata entries
    if nndb > 0:
        cdelays = "nodelay"
        for kdx, k in enumerate(lnd):
            for i in range(nrow):
                for j in range(ncol):
                    # skip constant head cells
                    if k == 0 and i == nrow - 1 and j in ccol:
                        continue
                    tag = f"{k + 1:02d}_{i + 1:02d}_{j + 1:02d}"
                    # create nodelay entry
                    # no delay beds
                    b = thicknd0[kdx]
                    d = [
                        ibcno,
                        (k, i, j),
                        cdelays,
                        hc[idx],
                        b,
                        1.0,
                        ccnd0[kdx],
                        crnd0[kdx],
                        theta,
                        999.0,
                        -999.0,
                        tag,
                    ]
                    sub6.append(d)
                    ibcno += 1

    # create delay bed packagedata entries and coarse-grained material storage
    S = []
    ndb = ndelaybeds[idx]
    if ndb > 0:
        cdelays = "delay"
        for kdx, k in enumerate(ldnd):
            for i in range(nrow):
                for j in range(ncol):
                    # skip constant head cells
                    if k == 0 and i == nrow - 1 and j in ccol:
                        continue
                    tag = f"{k + 1:02d}_{i + 1:02d}_{j + 1:02d}"
                    # create nodelay entry
                    d = [
                        ibcno,
                        (k, i, j),
                        cdelays,
                        dhc[kdx],
                        dz[kdx],
                        rnb[kdx],
                        cc,
                        cr,
                        theta,
                        kv,
                        dstart[kdx][i, j],
                        tag,
                    ]
                    sub6.append(d)
                    ibcno += 1

        # create S for aquifer, delay beds, and no-delay beds
        for k in range(nlay):
            sst = (1.0 - ffrac[k]) * zthick[k] * ss[k]
            S.append(sst)
    else:
        # create S for aquifer and no-delay beds
        for k in range(nlay):
            sst = (zthick[k] - thicknd0[k]) * ss[k]
            S.append(sst)

    maxibc = len(sub6)

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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

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
        icelltype=laytyp,
        cvoptions=cvopt[idx],
        k=hk,
        k33=vka,
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        # ss=S, sy=sy,
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
    # csub files
    opth = f"{name}.csub.obs"
    ibcsv = f"{name}.ib.strain.csv"
    skcsv = f"{name}.sk.strain.csv"
    copth = f"{name}.compaction.gridbin"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        print_input=True,
        boundnames=True,
        head_based=True,
        effective_stress_lag=True,
        specified_initial_interbed_state=True,
        save_flows=True,
        strainib_filerecord=ibcsv,
        straincg_filerecord=skcsv,
        compaction_filerecord=copth,
        ndelaycells=ndelaycells[idx],
        ninterbeds=maxibc,
        beta=0.0,
        cg_ske_cr=ss,
        packagedata=sub6,
    )
    orecarray = {}
    orecarray["csub_obs.csv"] = [
        ("tcomp1", "interbed-compaction", (44,)),
        ("tcomp2", "interbed-compaction", (140,)),
        ("tcomp3", "interbed-compaction", (240,)),
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
    sim = get_model(idx, test.workspace)

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
    fpth = os.path.join(test.workspace, cmppth, "csub_obs.csv")
    try:
        tc0 = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    diff = tc["TCOMP3"] - tc0["TCOMP3"]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-compaction difference ({diffmax}) "

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
    cbc_bud = ["CSUB-ELASTIC", "CSUB-INELASTIC"]
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
            for vv in v["q"]:
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
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        htol=htol[idx],
        compare="mf6_regression",
    )
    test.run()

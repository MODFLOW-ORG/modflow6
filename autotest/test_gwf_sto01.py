import os

import flopy
import numpy as np
import pytest
from conftest import try_get_target
from framework import TestFramework

cases = ["gwf_sto01"]
cmppth = "mfnwt"
tops = [0.0]
htol = [None for _ in range(len(cases))]
dtol = 1e-3
budtol = 1e-2
bud_lst = ["STO-SS_IN", "STO-SS_OUT", "STO-SY_IN", "STO-SY_OUT"]

# static model data
# temporal discretization
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
# tsmult = [1.0] + [1.0 for i in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# spatial discretization data
nlay, nrow, ncol = 3, 10, 10
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1000.0, 2000.0
botm = [-100, -150.0, -350.0]
strt = 0.0
hnoflo = 1e30
hdry = -1e30

# calculate hk
hk1fact = 1.0 / 50.0
hk1 = np.ones((nrow, ncol), dtype=float) * 0.5 * hk1fact
hk1[0, :] = 1000.0 * hk1fact
hk1[-1, :] = 1000.0 * hk1fact
hk1[:, 0] = 1000.0 * hk1fact
hk1[:, -1] = 1000.0 * hk1fact
hk = [20.0, hk1, 5.0]

# calculate vka
vka = [1e6, 7.5e-5, 1e6]

# all cells are active and layer 1 is convertible
ib = 1
laytyp = [1, 0, 0]

# solver options
nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0
newtonoptions = "NEWTON"
imsla = "BICGSTAB"

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

# storage and compaction data
ske = [6e-4, 3e-4, 6e-4]


# variant SUB package problem 3
def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    top = tops[idx]
    zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
    elevs = [top] + botm

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions=newtonoptions, save_flows=True
    )

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
        ss=ske,
        sy=0,
        storagecoefficient=None,
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

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    # build MODFLOW-NWT files
    cpth = cmppth
    ws = os.path.join(test.workspace, cpth)
    mc = flopy.modflow.Modflow(
        name,
        model_ws=ws,
        version=cpth,
        exe_name=try_get_target(test.targets, "mfnwt"),
    )
    dis = flopy.modflow.ModflowDis(
        mc,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        tsmult=tsmult,
        steady=steady,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    bas = flopy.modflow.ModflowBas(mc, ibound=ib, strt=strt, hnoflo=hnoflo, stoper=0.01)
    upw = flopy.modflow.ModflowUpw(
        mc,
        laytyp=laytyp,
        ipakcb=1001,
        hk=hk,
        vka=vka,
        ss=ske,
        sy=0.0,
        hdry=hdry,
    )
    chd = flopy.modflow.ModflowChd(mc, stress_period_data=cd)
    rch = flopy.modflow.ModflowRch(mc, rech=rech)
    wel = flopy.modflow.ModflowWel(mc, stress_period_data=wd)
    oc = flopy.modflow.ModflowOc(
        mc,
        stress_period_data=None,
        save_every=1,
        save_types=["save head", "save budget"],
    )
    fluxtol = (float(nlay * nrow * ncol) - 4.0) * rclose
    nwt = flopy.modflow.ModflowNwt(
        mc,
        headtol=hclose,
        fluxtol=fluxtol,
        maxiterout=nouter,
        linmeth=2,
        maxitinner=ninner,
        unitnumber=132,
        options="SPECIFIED",
        backflag=0,
        idroptol=0,
    )
    return sim, mc


def check_output(idx, test):
    # get results from listing file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.lst")
    budl = flopy.utils.Mf6ListBudget(fpth)
    names = list(bud_lst)
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]

    # get results from cbc file
    cbc_bud = ["STO-SS", "STO-SY"]
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


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        htol=htol[idx],
    )
    test.run()

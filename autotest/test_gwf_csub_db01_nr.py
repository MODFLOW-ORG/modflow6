import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = (
    "csub_db01a",
    "csub_db01b",
    "csub_db01c",
    "csub_db01d",
    "csub_db01e",
    "csub_db01f",
)

newtons = (True, True, True, True, True, True)
stress_lag = (None, True, None, True, None, True)
elastic = (True, True, False, False, False, False)
convertible = (True, True, True, True, False, False)
htol = None
dtol = 1e-3
budtol = 0.01
bud_lst = (
    "CSUB-CGELASTIC_IN",
    "CSUB-CGELASTIC_OUT",
)

# static model data
nlay, nrow, ncol = 2, 1, 2
nper = 3
tsp0 = 1.0
perlen = [tsp0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [200 for _ in range(nper - 1)]
tsmult = [1.0] + [1.0 for _ in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
delr, delc = 1000.0, 1000.0
top = 0.0
botm = [-10.0, -20.0]
zthick = [top - botm[0], botm[0] - botm[1]]
strt = -5.0
hmin = -15.0

shape3d = (nlay, nrow, ncol)

# npf variables
hk = 1e4
laytyp = 1
sy = 0.0

nouter, ninner = 50, 300
hclose, rclose, relax = 1e-9, 1e-3, 1.0

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# all cells are active
ib = 1

# chd data
ts_name = "CHD_TS"
c6 = [[1, 0, 0, ts_name]]
cd6 = {0: c6}
maxchd = len(cd6[0])

# static sub data
sgm = 1.7
sgs = 2.0
void = 0.82
theta = void / (1.0 + void)
beta = 0.0
crnd0 = 0.0
cc, cr = 6e-4, 6e-6
ini_stress = 0.0
kv = 1e-3
H0 = 0.0


def build_models(idx, test):
    name = cases[idx]
    newton = newtons[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    if newton:
        newtonoptions = "NEWTON"
        imsla = "BICGSTAB"
        rewet_record = None
        wetdry = None
    else:
        newtonoptions = None
        imsla = "CG"
        rewet_record = ["wetfct", 0.1, "iwetit", 1, "ihdwet", 0]
        wetdry = [1, 0]

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=f"{rclose} strict",
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
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

    obs_recarray = {
        "gwf_obs.csv": [
            ("h2_1_1", "HEAD", (1, 0, 0)),
            ("h1_1_2", "HEAD", (0, 0, 1)),
            ("h2_1_2", "HEAD", (1, 0, 1)),
        ]
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        gwf,
        pname="head_obs",
        digits=10,
        print_input=True,
        continuous=obs_recarray,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=laytyp,
        k=hk,
        rewet_record=rewet_record,
        wetdry=wetdry,
    )
    # storage
    if convertible[idx]:
        iconvert = 1
    else:
        iconvert = 0
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=iconvert,
        ss=0.0,
        sy=sy,
        steady_state={0: True},
        transient={1: True},
    )

    # create chd time series
    chnam = f"{name}.ch.ts"
    chd_ts = [
        (0.0, strt),
        (tsp0, strt),
        (365.25 + tsp0, hmin),
        (np.array(perlen).sum(), strt),
    ]

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
    )

    # initialize time series
    chd.ts.initialize(
        filename=chnam,
        timeseries=chd_ts,
        time_series_namerecord=[ts_name],
        interpolation_methodrecord=["linear"],
    )

    # csub observations
    obspos = [
        (0, 0, 1),
        (1, 0, 1),
    ]
    obstype = [
        "compaction-cell",
        "csub-cell",
    ]
    obstag = [
        "tcomp",
        "csub",
    ]
    obsarr = []
    for iobs, cobs in enumerate(obstype):
        for jobs, otup in enumerate(obspos):
            otag = f"{obstag[iobs]}{jobs + 1}"
            obsarr.append((otag, cobs, otup))
    obsloc = [
        (0, 0),
        (0, 9),
        (0, 18),
        (1, 0),
        (1, 9),
        (1, 18),
    ]
    obstag = [
        "dh1-1",
        "dh1-10",
        "dh1-19",
        "dh2-1",
        "dh2-10",
        "dh2-19",
    ]
    for loc, tag in zip(obsloc, obstag):
        obsarr.append((tag, "delay-head", loc))

    # csub files
    swt6 = []
    ibcno = 0
    if elastic[idx]:
        cct, crt = cr, cr
    else:
        cct, crt = cc, cr
    for k in range(nlay):
        d = [
            ibcno,
            (k, 0, 1),
            "delay",
            ini_stress,
            zthick[k],
            1.0,
            cct,
            crt,
            theta,
            kv,
            H0,
        ]
        swt6.append(d)
        ibcno += 1

    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        observations={"csub_obs.csv": obsarr},
        effective_stress_lag=stress_lag[idx],
        save_flows=True,
        beta=beta,
        sgm=sgm,
        sgs=sgs,
        ninterbeds=len(swt6),
        cg_theta=theta,
        cg_ske_cr=crnd0,
        packagedata=swt6,
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
    return sim, None


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
    cbc_bud = ("CSUB-CGELASTIC",)
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
        htol=htol,
    )
    test.run()

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

simname = "gwfcsubsk01"
cases = [f"{simname}a", f"{simname}b", f"{simname}c"]
dtol = 1e-3
budtol = 0.01
bud_lst = [
    "CSUB-CGELASTIC_IN",
    "CSUB-CGELASTIC_OUT",
    "CSUB-WATERCOMP_IN",
    "CSUB-WATERCOMP_OUT",
]
cvopt = [None, None, None]
constantcv = [True, True, True]
ndelaybeds = [0, 0, 0]
top = [0, 0, 15]
newton = [False, True, True]
htol = [None, None, 0.3]


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    cmp = get_model(idx, test.workspace / "mf6_regression")
    return sim, cmp


def get_model(idx, workspace):
    name = cases[idx]
    newtonoptions = None
    imsla = "CG"
    if newton[idx]:
        newtonoptions = "NEWTON"
        imsla = "BICGSTAB"

    # static model data
    nlay, nrow, ncol = 3, 10, 10
    nper = 31
    perlen = [1.0] + [365.2500000 for i in range(nper - 1)]
    nstp = [1] + [6 for i in range(nper - 1)]
    tsmult = [1.0] + [1.3 for i in range(nper - 1)]
    steady = [True] + [False for i in range(nper - 1)]
    delr, delc = 1000.0, 2000.0
    botm = [-100, -150.0, -350.0]
    zthick = [top[idx] - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
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
    sy = 0.0  # [0.1, 0., 0.]

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
    void = 0.82
    theta = void / (1.0 + void)
    sw = 4.65120000e-10 * 9806.65000000 * theta

    # no delay bed data
    nndb = 3
    lnd = [0, 1, 2]
    hc = [botm[-1] for k in range(nlay)]
    thicknd0 = [zthick[0], zthick[1], zthick[2]]
    ccnd0 = [6e-6, 3e-6, 6e-6]
    crnd0 = [6e-6, 3e-6, 6e-6]
    sfv = []
    sfe = []
    for k in range(nlay):
        sfv.append(ccnd0[k] * thicknd0[k])
        sfe.append(crnd0[k] * thicknd0[k])

    # sub output data
    ds15 = [0, 0, 0, 2052, 0, 0, 0, 0, 0, 0, 0, 0]
    ds16 = [0, nper - 1, 0, nstp[-1] - 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=str(workspace),
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
        top=top[idx],
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
        cvoptions=cvopt[idx],
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
    # csub files
    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        head_based=True,
        save_flows=True,
        ninterbeds=0,
        cg_theta=theta,
        cg_ske_cr=crnd0,
        packagedata=None,
    )
    obspos = [(0, 4, 4), (1, 4, 4), (2, 4, 4)]
    obstype = ["compaction-cell", "csub-cell"]
    obstag = ["tcomp", "csub"]
    obsarr = []
    for iobs, cobs in enumerate(obstype):
        for jobs, otup in enumerate(obspos):
            otag = f"{obstag[iobs]}{jobs + 1}"
            obsarr.append((otag, cobs, otup))

    obsarr2 = []
    obstype2 = [
        "csub",
        "inelastic-csub",
        "elastic-csub",
        "sk",
        "ske",
        "thickness",
        "theta",
        "interbed-compaction",
        "inelastic-compaction",
        "elastic-compaction",
        "delay-flowtop",
        "delay-flowbot",
    ]
    iobs = 0
    for cobs in obstype2:
        iobs += 1
        otag = f"obs{iobs:03d}"
        obsarr2.append((otag, cobs, (0,)))

    obstype3 = [
        "delay-preconstress",
        "delay-head",
        "delay-gstress",
        "delay-estress",
        "delay-compaction",
        "delay-thickness",
        "delay-theta",
    ]
    for cobs in obstype3:
        iobs += 1
        otag = f"obs{iobs:03d}"
        obsarr2.append((otag, cobs, (0,), (0,)))

    obsarr3 = []
    obstype4 = [
        "gstress-cell",
        "estress-cell",
        "thickness-cell",
        "coarse-csub",
        "wcomp-csub-cell",
        "coarse-compaction",
        "coarse-theta",
        "coarse-thickness",
        "csub-cell",
        "ske-cell",
        "sk-cell",
        "theta-cell",
        "compaction-cell",
    ]
    for cobs in obstype4:
        iobs += 1
        otag = f"obs{iobs:03d}"
        obsarr3.append((otag, cobs, obspos[-1]))

    orecarray = {}
    orecarray["csub_obs.csv"] = obsarr
    orecarray["interbed_obs.csv"] = obsarr2
    orecarray["coarse_cell_obs.csv"] = obsarr3

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


def check_output(idx, test):
    # MODFLOW 6 total compaction results
    fpth = os.path.join(test.workspace, "csub_obs.csv")
    tc = np.genfromtxt(fpth, names=True, delimiter=",")

    # regression compaction results
    cpth = "mf6_regression"
    fpth = os.path.join(test.workspace, cpth, "csub_obs.csv")
    tc0 = np.genfromtxt(fpth, names=True, delimiter=",")

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
                for j, key in enumerate(bud_lst):
                    line += f"{key + '_LST':>25s}"
                    line += f"{key + '_CBC':>25s}"
                    line += f"{key + '_DIF':>25s}"
                f.write(line + "\n")
            line = f"{d['totim'][i]:10g}"
            for j, key in enumerate(bud_lst):
                line += f"{d0[key][i]:25g}"
                line += f"{d[key][i]:25g}"
                line += f"{diff[i, j]:25g}"
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
        compare="mf6_regression",
    )
    test.run()

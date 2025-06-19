"""
This autotest is based on the MOC3D problem 1 autotest except that it
tests the zero-order decay for a simple one-dimensional flow problem.
The test ensures that concentrations do not go below zero (they do go
slightly negative but, it does ensure that the decay rate shuts off
where concentrations are zero.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "moc3d01zoda",
    "moc3d01zodb",
    "moc3d01zodc",
    "moc3d01zodd",
]
retardation = [None, 40, None, 40]
decay = [0.01, 0.01, 0.1, 0.1]
ist_package = [False, False, True, True]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 122, 1
    nper = 1
    perlen = [120]
    nstp = [240]
    tsmult = [1.0]
    steady = [True]
    delr = 0.1
    delc = 0.1
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 0.01
    laytyp = 0
    diffc = 0.0
    alphal = 0.1
    # ss = 0.
    # sy = 0.1

    nouter, ninner = 200, 300
    hclose, rclose, relax = 1e-8, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        # continue_=True,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
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
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        save_specific_discharge=True,
        icelltype=laytyp,
        k=hk,
        k33=hk,
    )
    # storage
    # sto = flopy.mf6.ModflowGwfsto(gwf, save_flows=False,
    #                              iconvert=laytyp[idx],
    #                              ss=ss[idx], sy=sy[idx],
    #                              steady_state={0: True, 2: True},
    #                              transient={1: True})

    # chd files
    c = {0: [[(0, 121, 0), 0.0000000]]}
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=c, save_flows=False, pname="CHD-1"
    )

    # wel files
    w = {0: [[(0, 0, 0), 0.001, 1.0]]}
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=w,
        save_flows=False,
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(
        sim,
        modelname=gwtname,
        save_flows=True,
        model_nam_file=f"{gwtname}.nam",
    )

    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(
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
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    strt = np.zeros((nlay, nrow, ncol))
    strt[0, 0, 0] = 0.0
    ic = flopy.mf6.ModflowGwtic(gwt, strt=strt, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="tvd", filename=f"{gwtname}.adv")

    # dispersion
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        diffc=diffc,
        alh=alphal,
        alv=alphal,
        ath1=0.0,
        atv=0.0,
        filename=f"{gwtname}.dsp",
    )

    # storage
    theta_mobile = 0.1  # vol mobile voids per cell volume
    volfrac_immobile = 0.0
    theta_immobile = 0.0
    if ist_package[idx]:
        # if dual domain, then assume half of cell is mobile and other half is immobile
        volfrac_immobile = 0.5
        theta_immobile = theta_mobile
        porosity_immobile = theta_immobile / volfrac_immobile
    volfrac_mobile = 1.0 - volfrac_immobile
    porosity_mobile = theta_mobile / volfrac_mobile

    rtd = retardation[idx]
    sorption = None
    kd = None
    rhob = None
    if rtd is not None:
        rhob = 1.0
        kd = (rtd - 1.0) * theta_mobile / rhob
        rhobm = rhob
        sorption = "linear"

    decay_rate = decay[idx]
    zero_order_decay = False
    if decay_rate is not None:
        zero_order_decay = True

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity_mobile,
        zero_order_decay=zero_order_decay,
        decay=decay_rate,
        decay_sorbed=decay_rate,
        sorption=sorption,
        distcoef=kd,
        bulk_density=rhob,
    )

    if ist_package[idx]:
        ist = flopy.mf6.ModflowGwtist(
            gwt,
            cim_filerecord=f"{gwtname}.ist.ucn",
            sorption=sorption,
            zero_order_decay=True,
            cim=0.0,
            volfrac=volfrac_immobile,
            porosity=porosity_immobile,
            zetaim=1.0,
            decay=decay_rate,
            bulk_density=rhob,
            distcoef=kd,
            decay_sorbed=decay_rate,
        )

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def make_plot_ct(tssim, fname=None):
    """Concentration versus time plot"""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1)
    mec = ["red", "blue", "green"]
    iskip = 2
    tssim = tssim[::iskip]
    for i, l in enumerate(["x=0.05", "x=4.05", "x=11.05"]):
        ax.plot(
            tssim[:, 0],
            tssim[:, i + 1],
            marker="o",
            ls="none",
            mec=mec[i],
            mfc="none",
            markersize="4",
            label=l,
        )

    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Normalized Concentration, dimensionless")
    plt.legend()

    if fname is not None:
        plt.savefig(fname, bbox_inches="tight")
    return


def make_plot_cd(cobj, fname=None):
    """Concentration versus distance plot"""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1)
    delr = 0.1
    system_length = 12.0
    ncol = 122
    iskip = 1
    mec = ["red", "blue", "green"]
    x = np.linspace(0.5 * delr, system_length - 0.5 * delr, ncol)
    for i, t in enumerate([6.0, 60.0, 120.0]):
        conc = cobj.get_data(totim=t).flatten()
        ax.plot(
            x[::iskip],
            conc[::iskip],
            marker="o",
            ls="none",
            mec=mec[i],
            mfc="none",
            markersize="4",
            label=f"t={t} s",
        )

    ax.set_xlabel("Distance (cm)")
    ax.set_ylabel("Normalized Concentration, dimensionless")
    plt.legend()

    if fname is not None:
        plt.savefig(fname, bbox_inches="tight")
    return


def plot_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name
    sim = test.sims[0]
    gwt = sim.get_model(gwtname)
    cobj = gwt.output.concentration()
    station = [(0, 0, 0), (0, 40, 0), (0, 110, 0)]
    tssim = cobj.get_ts(station)

    # concentration versus time
    fname = test.workspace / "fig-ct.pdf"
    make_plot_ct(tssim, fname)

    # concentration versus distance
    fname = test.workspace / "fig-cd.pdf"
    make_plot_cd(cobj, fname)


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name
    sim = test.sims[0]
    gwt = sim.get_model(gwtname)
    cobj = gwt.output.concentration()
    bobj = gwt.output.budget()
    station = [(0, 0, 0), (0, 40, 0), (0, 110, 0)]
    tssim = cobj.get_ts(station)

    # Check to make sure decay rates in budget file are correct.  If there is
    # enough mass in the cell, then the qdecay value in the budget file
    # should be equal to decay_rate * vcell * porosity
    records = bobj.get_data(text="decay")
    qdecay_budfile = records[0].flatten()
    conc = cobj.get_data().flatten()
    delt = 0.5
    vcell = 0.1 * 0.1 * 1.0
    porosity = 0.1
    decay_rate = decay[idx]
    for i in range(122):
        if conc[i] / delt > decay_rate:
            qknown = -decay_rate * vcell * porosity
            errmsg = (
                f"Decay rate in budget file for cell {i} should be "
                f"{qdecay_budfile[i]} but found {qknown} instead."
            )
            assert np.allclose(qdecay_budfile[i], qknown), errmsg
        # print(i, qdecay_budfile[i], conc[i])

    # get immobile domain concentration object
    fpth = os.path.join(test.workspace, f"{gwtname}.ist.ucn")
    cimobj = None
    if os.path.isfile(fpth):
        try:
            cimobj = flopy.utils.HeadFile(fpth, precision="double", text="CIM")
        except:
            assert False, f'could not load data from "{fpth}"'

        records = bobj.get_data(text="immobile domain")
        qim_budfile = records[0]["q"]
        cim = cimobj.get_data().flatten()
        cm = cobj.get_data().flatten()
        zetaim = 1.0
        qim_calculated = (cim - cm) * zetaim * vcell
        # for i in range(122):
        #    print(i, cm[i], cim[i], qim_budfile[i], qim_calculated[i])
        errmsg = (
            "Mass transfer rates from the gwt budget file do not "
            "compare with mass transfer rates calculated from "
            "simulated mobile and immobile domain concentrations\n"
            f"{qim_budfile} /= {qim_calculated}"
        )
        np.allclose(qim_budfile, qim_calculated), errmsg

    # compare every tenth time
    tssim = tssim[::10]
    # print(tssim)

    # answer for case with decay and no sorption; taken from run that appeared
    # to have the correct answer.
    tsresa = [
        [5.00000000e-01, 2.99857201e-01, -2.56702053e-18, -6.21538564e-43],
        [5.50000000e00, 9.27984664e-01, -9.48839600e-13, -3.68700082e-33],
        [1.05000000e01, 9.75342997e-01, -3.14583685e-10, -1.94112207e-27],
        [1.55000000e01, 9.82978097e-01, -6.67553006e-09, -2.47887998e-23],
        [2.05000000e01, 9.84532610e-01, -4.94300402e-08, -3.45215269e-20],
        [2.55000000e01, 9.84885335e-01, -1.43321212e-07, -1.02381037e-17],
        [3.05000000e01, 9.84970692e-01, 8.74919559e-08, -9.62720397e-16],
        [3.55000000e01, 9.84992281e-01, 3.29076093e-02, -3.78602524e-14],
        [4.05000000e01, 9.84997920e-01, 1.57687461e-01, -7.68130469e-13],
        [4.55000000e01, 9.84999429e-01, 3.01668660e-01, -9.37217658e-12],
        [5.05000000e01, 9.84999841e-01, 4.18583229e-01, -7.55021768e-11],
        [5.55000000e01, 9.84999955e-01, 4.96262647e-01, -4.24745302e-10],
        [6.05000000e01, 9.84999987e-01, 5.41229549e-01, -1.74881406e-09],
        [6.55000000e01, 9.84999996e-01, 5.64736374e-01, -5.50775463e-09],
        [7.05000000e01, 9.84999999e-01, 5.76095202e-01, -1.37606603e-08],
        [7.55000000e01, 9.85000000e-01, 5.81252028e-01, -2.81572476e-08],
        [8.05000000e01, 9.85000000e-01, 5.83478177e-01, -4.88424746e-08],
        [8.55000000e01, 9.85000000e-01, 5.84400410e-01, -7.48093720e-08],
        [9.05000000e01, 9.85000000e-01, 5.84769691e-01, -1.06424623e-07],
        [9.55000000e01, 9.85000000e-01, 5.84913425e-01, -1.42793519e-07],
        [1.00500000e02, 9.85000000e-01, 5.84968051e-01, -1.77377488e-07],
        [1.05500000e02, 9.85000000e-01, 5.84988395e-01, -2.05441285e-07],
        [1.10500000e02, 9.85000000e-01, 5.84995843e-01, -2.25227219e-07],
        [1.15500000e02, 9.85000000e-01, 5.84998528e-01, -2.34638026e-07],
    ]

    # answer for case with decay and sorption
    tsresb = [
        [5.00000000e-01, 1.09202432e-02, -1.40159170e-68, -2.09981357e-181],
        [5.50000000e00, 1.09676704e-01, -9.20100316e-60, -1.63693421e-168],
        [1.05000000e01, 1.92151205e-01, -1.19261506e-54, -8.57633085e-160],
        [1.55000000e01, 2.61771923e-01, -6.41130367e-51, -7.48882965e-153],
        [2.05000000e01, 3.21147740e-01, -5.74798815e-48, -5.20924225e-147],
        [2.55000000e01, 3.72249517e-01, -1.60821228e-45, -6.14790781e-142],
        [3.05000000e01, 4.16576496e-01, -1.97762152e-43, -1.91416424e-137],
        [3.55000000e01, 4.55290693e-01, -1.31524619e-41, -2.09339517e-133],
        [4.05000000e01, 4.89308401e-01, -5.41319756e-40, -9.80957395e-130],
        [4.55000000e01, 5.19362011e-01, -1.51258771e-38, -2.27640662e-126],
        [5.05000000e01, 5.46043609e-01, -3.06808692e-37, -2.91831182e-123],
        [5.55000000e01, 5.69836361e-01, -4.75075541e-36, -2.24992540e-120],
        [6.05000000e01, 5.91137653e-01, -5.83989930e-35, -1.11593389e-117],
        [6.55000000e01, 6.10276851e-01, -5.87961560e-34, -3.76028096e-115],
        [7.05000000e01, 6.27529340e-01, -4.97294419e-33, -9.00232238e-113],
        [7.55000000e01, 6.43127364e-01, -3.60829888e-32, -1.58932043e-110],
        [8.05000000e01, 6.57268283e-01, -2.28577909e-31, -2.13497151e-108],
        [8.55000000e01, 6.70120926e-01, -1.28308590e-30, -2.24104025e-106],
        [9.05000000e01, 6.81830520e-01, -6.46359215e-30, -1.88054775e-104],
        [9.55000000e01, 6.92522567e-01, -2.95413799e-29, -1.28661098e-102],
        [1.00500000e02, 7.02305955e-01, -1.23662165e-28, -7.30091788e-101],
        [1.05500000e02, 7.11275477e-01, -4.78052344e-28, -3.48810622e-99],
        [1.10500000e02, 7.19513925e-01, -1.71903432e-27, -1.42175219e-97],
        [1.15500000e02, 7.27093862e-01, -5.78658515e-27, -5.00225543e-96],
    ]

    # answer for case with decay and immobile decay
    tsresc = [
        [5.00000000e-01, 1.77197158e-01, -4.98160896e-26, -3.68736386e-61],
        [5.50000000e00, 6.13562680e-01, -1.95022911e-18, -1.41020506e-49],
        [1.05000000e01, 6.79902159e-01, -5.89336159e-15, -9.07716595e-43],
        [1.55000000e01, 6.95851560e-01, -8.67651988e-13, -9.85058240e-38],
        [2.05000000e01, 7.00511382e-01, -2.61536943e-11, -9.44216128e-34],
        [2.55000000e01, 7.01967536e-01, -2.99804604e-10, -1.72760736e-30],
        [3.05000000e01, 7.02428631e-01, -1.81199372e-09, -9.37792675e-28],
        [3.55000000e01, 7.02574913e-01, -7.01024063e-09, -2.00235899e-25],
        [4.05000000e01, 7.02621325e-01, -1.96671395e-08, -2.04073887e-23],
        [4.55000000e01, 7.02636049e-01, -4.35207683e-08, -1.14127936e-21],
        [5.05000000e01, 7.02640721e-01, -8.04787760e-08, -3.88681552e-20],
        [5.55000000e01, 7.02642203e-01, -1.29518117e-07, -8.73292857e-19],
        [6.05000000e01, 7.02642673e-01, -1.86905909e-07, -1.37854507e-17],
        [6.55000000e01, 7.02642822e-01, -2.47512354e-07, -1.60793164e-16],
        [7.05000000e01, 7.02642869e-01, -3.06381734e-07, -1.44380769e-15],
        [7.55000000e01, 7.02642884e-01, -3.59816171e-07, -1.03235722e-14],
        [8.05000000e01, 7.02642889e-01, -4.05737318e-07, -6.04609190e-14],
        [8.55000000e01, 7.02642890e-01, -4.43503801e-07, -2.96999301e-13],
        [9.05000000e01, 7.02642891e-01, -4.73484022e-07, -1.24861521e-12],
        [9.55000000e01, 7.02642891e-01, -4.96613642e-07, -4.57071104e-12],
        [1.00500000e02, 7.02642891e-01, -5.14050080e-07, -1.47865421e-11],
        [1.05500000e02, 7.02642891e-01, -5.26949080e-07, -4.28208988e-11],
        [1.10500000e02, 7.02642891e-01, -5.36344676e-07, -1.12254606e-10],
        [1.15500000e02, 7.02642891e-01, -5.43101091e-07, -2.69001197e-10],
    ]

    # answer for case with decay, sorption and immobile decay
    tsresd = [
        [5.00000000e-01, 7.91464499e-03, -7.94470260e-63, -9.01303292e-162],
        [5.50000000e00, 3.77321075e-02, -8.12858000e-54, -9.48147419e-149],
        [1.05000000e01, 4.14136035e-02, -4.43077646e-49, -1.82848874e-140],
        [1.55000000e01, 4.18681499e-02, -9.43913042e-46, -5.02351377e-134],
        [2.05000000e01, 4.19242717e-02, -3.50882796e-43, -1.06515915e-128],
        [2.55000000e01, 4.19312010e-02, -4.33956148e-41, -3.86156404e-124],
        [3.05000000e01, 4.19320565e-02, -2.51511261e-39, -3.78368013e-120],
        [3.55000000e01, 4.19321621e-02, -8.37857960e-38, -1.34435242e-116],
        [4.05000000e01, 4.19321752e-02, -1.82963426e-36, -2.12027749e-113],
        [4.55000000e01, 4.19321768e-02, -2.86419101e-35, -1.71793239e-110],
        [5.05000000e01, 4.19321770e-02, -3.42542504e-34, -7.97777855e-108],
        [5.55000000e01, 4.19321770e-02, -3.27974550e-33, -2.31012283e-105],
        [6.05000000e01, 4.19321770e-02, -2.60484546e-32, -4.45773639e-103],
        [6.55000000e01, 4.19321770e-02, -1.76395035e-31, -6.04610029e-101],
        [7.05000000e01, 4.19321770e-02, -1.04089801e-30, -6.01997960e-99],
        [7.55000000e01, 4.19321770e-02, -5.44705062e-30, -4.56109021e-97],
        [8.05000000e01, 4.19321770e-02, -2.56431053e-29, -2.70975328e-95],
        [8.55000000e01, 4.19321770e-02, -1.09898545e-28, -1.29473714e-93],
        [9.05000000e01, 4.19321770e-02, -4.33058298e-28, -5.08390429e-92],
        [9.55000000e01, 4.19321770e-02, -1.58231318e-27, -1.67119010e-90],
        [1.00500000e02, 4.19321770e-02, -5.39953641e-27, -4.67335620e-89],
        [1.05500000e02, 4.19321770e-02, -1.73153803e-26, -1.12737319e-87],
        [1.10500000e02, 4.19321770e-02, -5.24634960e-26, -2.37493203e-86],
        [1.15500000e02, 4.19321770e-02, -1.50895159e-25, -4.41626189e-85],
    ]

    tsresa = np.array(tsresa)
    tsresb = np.array(tsresb)
    tsresc = np.array(tsresc)
    tsresd = np.array(tsresd)
    tsreslist = [tsresa, tsresb, tsresc, tsresd]
    tsres = tsreslist[idx]
    errmsg = (
        "Simulated concentrations do not match with known solution.\n"
        f"{tssim} /= {tsres}"
    )
    if tsres is not None:
        assert np.allclose(tsres, tssim), errmsg


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
    )
    test.run()

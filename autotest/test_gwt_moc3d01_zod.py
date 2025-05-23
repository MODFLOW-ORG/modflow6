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
        [5.00000000e-01, 2.80021420e-01, -6.21654981e-15, 2.94512187e-39],
        [5.50000000e+00, 9.02564301e-01, -4.84253930e-13, -1.09618214e-29],
        [1.05000000e+01, 9.67062210e-01, -1.05978896e-10, -6.67436392e-25],
        [1.55000000e+01, 9.78915405e-01, -3.53821341e-09, -8.36601978e-22],
        [2.05000000e+01, 9.81493335e-01, -2.99478585e-08, -1.17302104e-19],
        [2.55000000e+01, 9.82102870e-01, -9.51070742e-08, -4.59220902e-18],
        [3.05000000e+01, 9.82254646e-01, 8.60641012e-08, -1.00270370e-16],
        [3.55000000e+01, 9.82293842e-01, 3.30654200e-02, -1.87066152e-15],
        [4.05000000e+01, 9.82304244e-01, 1.57628457e-01, -2.73744064e-14],
        [4.55000000e+01, 9.82307064e-01, 3.01251363e-01, -4.59487047e-13],
        [5.05000000e+01, 9.82307841e-01, 4.17994177e-01, -6.37033069e-12],
        [5.55000000e+01, 9.82308058e-01, 4.95717338e-01, -5.73880104e-11],
        [6.05000000e+01, 9.82308119e-01, 5.40830243e-01, -3.51018542e-10],
        [6.55000000e+01, 9.82308137e-01, 5.64486282e-01, -1.51583463e-09],
        [7.05000000e+01, 9.82308142e-01, 5.75955359e-01, -4.83651532e-09],
        [7.55000000e+01, 9.82308143e-01, 5.81180311e-01, -1.19745875e-08],
        [8.05000000e+01, 9.82308144e-01, 5.83443818e-01, -2.40842909e-08],
        [8.55000000e+01, 9.82308144e-01, 5.84384825e-01, -4.09016940e-08],
        [9.05000000e+01, 9.82308144e-01, 5.84762931e-01, -6.04284362e-08],
        [9.55000000e+01, 9.82308144e-01, 5.84910599e-01, -7.98573160e-08],
        [1.00500000e+02, 9.82308144e-01, 5.84966906e-01, -9.94503185e-08],
        [1.05500000e+02, 9.82308144e-01, 5.84987944e-01, -1.22126257e-07],
        [1.10500000e+02, 9.82308144e-01, 5.84995668e-01, -1.42879932e-07],
        [1.15500000e+02, 9.82308144e-01, 5.84998462e-01, -1.59348081e-07],
    ]

    # answer for case with decay and sorption
    tsresb = [
        [5.00000000e-001, 1.08536585e-002, 1.50230289e-065, 7.03720789e-179],
        [5.50000000e000, 1.05972468e-001, -1.20770394e-056, -1.16040412e-164],
        [1.05000000e001, 1.81719287e-001, -1.73997465e-052, -2.69101712e-156],
        [1.55000000e001, 2.43729942e-001, -1.86616329e-049, -1.16363650e-149],
        [2.05000000e001, 2.95729413e-001, -7.29308335e-047, -3.99240889e-144],
        [2.55000000e001, 3.40178152e-001, -1.28823958e-044, -2.35667043e-139],
        [3.05000000e001, 3.78744937e-001, -1.07712892e-042, -3.79220669e-135],
        [3.55000000e001, 4.12607082e-001, -4.92494757e-041, -2.24709157e-131],
        [4.05000000e001, 4.42623555e-001, -1.40845582e-039, -6.05974088e-128],
        [4.55000000e001, 4.69439605e-001, -2.78724557e-038, -8.68988128e-125],
        [5.05000000e001, 4.93551649e-001, -4.10767481e-037, -7.43207001e-122],
        [5.55000000e001, 5.15348967e-001, -4.75649699e-036, -4.11388640e-119],
        [6.05000000e001, 5.35142116e-001, -4.50278413e-035, -1.55946013e-116],
        [6.55000000e001, 5.53183571e-001, -3.59069429e-034, -4.21269863e-114],
        [7.05000000e001, 5.69682313e-001, -2.46839024e-033, -8.36120511e-112],
        [7.55000000e001, 5.84813940e-001, -1.48978956e-032, -1.25121687e-109],
        [8.05000000e001, 5.98727767e-001, -8.01195046e-032, -1.44498742e-107],
        [8.55000000e001, 6.11551927e-001, -3.88652458e-031, -1.31567108e-105],
        [9.05000000e001, 6.23397138e-001, -1.71810925e-030, -9.63155190e-104],
        [9.55000000e001, 6.34359576e-001, -6.98234537e-030, -5.77117544e-102],
        [1.00500000e002, 6.44523151e-001, -2.62836345e-029, -2.87651068e-100],
        [1.05500000e002, 6.53961371e-001, -9.22475734e-029, -1.21009541e-098],
        [1.10500000e002, 6.62738943e-001, -3.03610362e-028, -4.35314284e-097],
        [1.15500000e002, 6.70913131e-001, -9.41863049e-028, -1.35492374e-095],
    ]

    # answer for case with decay and immobile decay
    tsresc = [
        [5.00000000e-01, 1.65743548e-01, -4.90064300e-23, -3.62743323e-58],
        [5.50000000e+00, 5.51785563e-01, -6.12543069e-19, -8.77247834e-48],
        [1.05000000e+01, 6.21289408e-01, -2.00557329e-16, -5.37346011e-42],
        [1.55000000e+01, 6.40753522e-01, -1.47570823e-14, -8.46420185e-38],
        [2.05000000e+01, 6.47093855e-01, -8.66843223e-13, -1.87670194e-34],
        [2.55000000e+01, 6.49242740e-01, -2.19755527e-11, -1.26768115e-31],
        [3.05000000e+01, 6.49974496e-01, -2.51505556e-10, -3.66835643e-29],
        [3.55000000e+01, 6.50223670e-01, -1.60029963e-09, -5.27338566e-27],
        [4.05000000e+01, 6.50308495e-01, -6.57525906e-09, -4.15193388e-25],
        [4.55000000e+01, 6.50337367e-01, -1.94556062e-08, -1.94561577e-23],
        [5.05000000e+01, 6.50347194e-01, -4.48948692e-08, -5.84308335e-22],
        [5.55000000e+01, 6.50350539e-01, -8.55687624e-08, -1.20138126e-20],
        [6.05000000e+01, 6.50351677e-01, -1.40383011e-07, -1.78545751e-19],
        [6.55000000e+01, 6.50352064e-01, -2.04424667e-07, -2.00166443e-18],
        [7.05000000e+01, 6.50352196e-01, -2.70732861e-07, -1.75198960e-17],
        [7.55000000e+01, 6.50352241e-01, -3.32744054e-07, -1.23148798e-16],
        [8.05000000e+01, 6.50352257e-01, -3.86026011e-07, -7.12703181e-16],
        [8.55000000e+01, 6.50352262e-01, -4.28731200e-07, -3.47140479e-15],
        [9.05000000e+01, 6.50352263e-01, -4.61082450e-07, -1.45946806e-14],
        [9.55000000e+01, 6.50352264e-01, -4.84506720e-07, -5.43265199e-14],
        [1.00500000e+02, 6.50352264e-01, -5.00868717e-07, -2.63522538e-13],
        [1.05500000e+02, 6.50352264e-01, -5.11977756e-07, -1.18779616e-12],
        [1.10500000e+02, 6.50352264e-01, -5.19353201e-07, -4.63830492e-12],
        [1.15500000e+02, 6.50352264e-01, -5.24163707e-07, -1.59516004e-11],
    ]

    # answer for case with decay, sorption and immobile decay
    tsresd = [
        [5.00000000e-001, 7.83710928e-003, -5.80158910e-057, -2.25314132e-155],
        [5.50000000e000, 3.62624591e-002, -6.24807157e-049, -4.24647039e-143],
        [1.05000000e001, 3.94430245e-002, -4.21132371e-045, -3.49914067e-135],
        [1.55000000e001, 3.97989039e-002, -1.88015521e-042, -4.22140911e-129],
        [2.05000000e001, 3.98387239e-002, -2.08000273e-040, -4.21326343e-124],
        [2.55000000e001, 3.98431794e-002, -9.73158626e-039, -7.70768515e-120],
        [3.05000000e001, 3.98436780e-002, -2.51979242e-037, -4.05412972e-116],
        [3.55000000e001, 3.98437337e-002, -4.22148524e-036, -8.16181323e-113],
        [4.05000000e001, 3.98437400e-002, -5.06137417e-035, -7.64733456e-110],
        [4.55000000e001, 3.98437407e-002, -4.65644754e-034, -3.83768947e-107],
        [5.05000000e001, 3.98437408e-002, -3.45690053e-033, -1.14526800e-104],
        [5.55000000e001, 3.98437408e-002, -2.14998235e-032, -2.20207186e-102],
        [6.05000000e001, 3.98437408e-002, -1.15268865e-031, -2.90478786e-100],
        [6.55000000e001, 3.98437408e-002, -5.44749494e-031, -2.76391812e-098],
        [7.05000000e001, 3.98437408e-002, -2.30981090e-030, -1.97569761e-096],
        [7.55000000e001, 3.98437408e-002, -8.91352546e-030, -1.09703810e-094],
        [8.05000000e001, 3.98437408e-002, -3.16720540e-029, -4.86540314e-093],
        [8.55000000e001, 3.98437408e-002, -1.04623376e-028, -1.76436266e-091],
        [9.05000000e001, 3.98437408e-002, -3.23874014e-028, -5.33676335e-090],
        [9.55000000e001, 3.98437408e-002, -9.45855357e-028, -1.36963650e-088],
        [1.00500000e002, 3.98437408e-002, -2.62075822e-027, -3.02677070e-087],
        [1.05500000e002, 3.98437408e-002, -6.92256308e-027, -5.83417810e-086],
        [1.10500000e002, 3.98437408e-002, -1.75036545e-026, -9.91958533e-085],
        [1.15500000e002, 3.98437408e-002, -4.25158985e-026, -1.50254693e-083],
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

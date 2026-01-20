"""
Test the IST Package with a one-dimensional flow problem
with dual-domain transport.  Flow is from left to right with
constant concentration equal to 1.0 for first 20 days and then
set to 0.0 for next 30 days.  The results are compared to
the results of an MT3D simulation sent by Sorab Panday.  The MT3D
results had many transport time steps, but they were interpolated
onto an even 1-day interval.  The mf6 results are also
interpolated onto the 1-day interval for comparison.  The test
passes if the difference in simulated concentration in column 300
between mf6 and mt3d is less than 0.05.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ist02"]
nlay, nrow, ncol = 1, 1, 300

mt3d_times = np.arange(1.0, 51.0, 1.0)
mt3d_conc = np.array(
    [
        0.000000e00,
        2.060000e-25,
        7.220000e-10,
        9.541440e-04,
        8.276250e-02,
        2.328700e-01,
        3.536184e-01,
        4.644795e-01,
        5.631907e-01,
        6.484800e-01,
        7.204300e-01,
        7.799700e-01,
        8.284300e-01,
        8.673400e-01,
        8.982000e-01,
        9.224300e-01,
        9.412600e-01,
        9.557800e-01,
        9.668900e-01,
        9.753300e-01,
        9.817000e-01,
        9.864900e-01,
        9.900600e-01,
        9.917600e-01,
        9.119200e-01,
        7.632500e-01,
        6.435600e-01,
        5.334700e-01,
        4.353300e-01,
        3.504500e-01,
        2.788000e-01,
        2.194700e-01,
        1.711600e-01,
        1.323600e-01,
        1.015800e-01,
        7.741230e-02,
        5.860264e-02,
        4.411496e-02,
        3.302846e-02,
        2.460321e-02,
        1.824100e-02,
        1.346449e-02,
        9.894626e-03,
        7.245270e-03,
        5.285890e-03,
        3.843070e-03,
        2.784940e-03,
        2.011888e-03,
        1.448659e-03,
        1.040850e-03,
    ]
)


def build_models(idx, test):
    perlen = [20.0, 30.0]
    nper = len(perlen)
    nstp = [100, 100]
    tsmult = [1.0, 1.0]
    delr = 0.5
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 9.5
    hk = 1000.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 0.97

    tdis_rc = []
    for id in range(nper):
        tdis_rc.append((perlen[id], nstp[id], tsmult[id]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    newtonoptions = None
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        newtonoptions=newtonoptions,
    )

    # create iterative model solution and register the gwf model with it
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=0, k=hk, k33=hk)

    # chd files
    chdspdict = {
        0: [[(0, 0, 0), 10.0], [(0, 0, ncol - 1), 9.0]],
    }
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        print_input=True,
        print_flows=True,
        stress_period_data=chdspdict,
        save_flows=False,
        pname="CHD-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)

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
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")

    # dispersion
    xt3d_off = True
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=xt3d_off,
        diffc=0.0,
        alh=0.5,
        ath1=0.0,
    )

    # mass storage and transfer
    bulk_density = 1.6
    thetam = 0.14
    thetaim = 0.12
    volfracm = thetam / (thetam + thetaim)
    volfracim = 1.0 - volfracm
    porositym = thetam / volfracm
    porosityim = thetaim / volfracim
    mst = flopy.mf6.ModflowGwtmst(
        gwt,
        sorption="LINEAR",
        porosity=porositym,
        bulk_density=bulk_density,
        distcoef=0.1,
    )

    # immobile storage and transfer
    cim_filerecord = f"{gwtname}.ist.ucn"
    ist = flopy.mf6.ModflowGwtist(
        gwt,
        sorption="LINEAR",
        save_flows=True,
        cim_filerecord=cim_filerecord,
        cim=0.0,
        porosity=porosityim,
        volfrac=volfracim,
        bulk_density=bulk_density,
        zetaim=0.1,
        distcoef=0.1,
    )

    # cnc
    cncspdict = {
        0: [[(0, 0, 0), 1.0]],
        1: [[(0, 0, 0), 0.0]],
    }
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        print_input=True,
        print_flows=True,
        maxbound=1,
        stress_period_data=cncspdict,
        save_flows=False,
        pname="CNC-1",
    )

    # sources
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    obs_data = {
        f"{gwtname}.obs.csv": [
            ("(1-1-300)", "CONCENTRATION", (0, 0, ncol - 1)),
        ],
    }
    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        # pname="conc_obs",
        filename=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
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


def plot_output(idx, test):
    name = test.name
    ws = test.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    output = gwt.obs.output
    obs_names = output.obs_names
    data = output.obs(f=obs_names[0]).data

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data["totim"], data["(1-1-300)"], "k-", label="mf6")
    ax.plot(mt3d_times, mt3d_conc, "ko", label="mt3d")
    plt.xlabel("time, in days")
    plt.ylabel("concentration, dimensionless")
    plt.legend()
    fname = os.path.join(ws, gwtname + ".png")
    plt.savefig(fname)


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name
    gwfname = "gwf_" + name

    # load the observed concentrations in column 300
    fname = os.path.join(test.workspace, gwtname + ".obs.csv")
    assert os.path.isfile(fname), f"file not found: {fname}"
    simvals = np.genfromtxt(fname, names=True, delimiter=",", deletechars="")

    # interpolate mf6 results to same times as mt3d
    mf6conc_interp = np.interp(mt3d_times, simvals["time"], simvals["(1-1-300)"])

    # calculate difference between mf6 and mt3d
    atol = 0.05
    diff = mf6conc_interp - mt3d_conc
    success = True
    print("index mf6 mt3d diff")
    for i in range(mf6conc_interp.shape[0]):
        print(f"{i} {mf6conc_interp[i]:.3f} {mt3d_conc[i]:.3f} {diff[i]:.3f}")
        if abs(diff[i]) > atol:
            success = False
    assert success, "Conc difference between mf6 and mt3d > 0.05"


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

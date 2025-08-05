"""
Test the immobile domain isotherms (Freundlich and Langmuir) in the IST Package
using a one-dimensional flow problem.  Compare the simulated sorbate
concentrations output from the model with sorbate concentrations calculated
in this script using python.  Problem is adapted from test_gwt_ist02.py.
"""

import os
import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ist03a", "ist03b", "ist03c"]
sorption_idx = ["LINEAR", "FREUNDLICH", "LANGMUIR"]
distcoef_idx = [0.1, 0.1, 0.1]
sp2_idx = [None, 0.7, 0.7]

nlay, nrow, ncol = 1, 1, 300


def build_models(idx, test):
    sorption = sorption_idx[idx]
    distcoef = distcoef_idx[idx]
    sp2 = sp2_idx[idx]

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
        sorbate_filerecord=f"{gwtname}.mst.csrb",
        sorption=sorption,
        porosity=porositym,
        bulk_density=bulk_density,
        distcoef=distcoef,
        sp2=sp2,
    )

    # immobile storage and transfer
    ist = flopy.mf6.ModflowGwtist(
        gwt,
        sorption=sorption,
        save_flows=True,
        cim_filerecord=f"{gwtname}.ist.ucn",
        sorbate_filerecord=f"{gwtname}.ist.csrb",
        cimprintrecord=[("COLUMNS", 5, "WIDTH", 11, "DIGITS", 4, "GENERAL")],
        cim=0.0,
        porosity=porosityim,
        volfrac=volfracim,
        bulk_density=bulk_density,
        zetaim=0.1,
        distcoef=distcoef,
        sp2=sp2,
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
    plot_title = sorption_idx[idx]
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

    # mobile aqueous sorbed concentrations
    fpth = pl.Path(ws) / f"{gwtname}.mst.csrb"
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="SORBATE")
    csrb = cobj.get_alldata().reshape(200, 300)

    # immobile aqueous
    fpth = pl.Path(ws) / f"{gwtname}.ist.ucn"
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CIM")
    cim = cobj.get_alldata().reshape(200, 300)

    # immobile sorbed
    fpth = pl.Path(ws) / f"{gwtname}.ist.csrb"
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="SORBATE")
    cimsrb = cobj.get_alldata().reshape(200, 300)

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(data["totim"], data["(1-1-300)"], "k-", label="mobile aqueous")
    ax.plot(data["totim"], csrb[:, 299], "b-", label="mobile sorbate")
    ax.plot(data["totim"], cim[:, 299], "k--", label="immobile aqueous")
    ax.plot(data["totim"], cimsrb[:, 299], "b--", label="immobile sorbate")
    plt.xlabel("time, in days")
    plt.ylabel("concentration")
    plt.legend()
    plt.title(plot_title)
    fname = os.path.join(ws, gwtname + ".png")
    plt.savefig(fname)


def check_output(idx, test):
    sorption = sorption_idx[idx]
    distcoef = distcoef_idx[idx]
    sp2 = sp2_idx[idx]

    name = test.name
    gwtname = "gwt_" + name
    gwfname = "gwf_" + name

    # load the observed concentrations in column 300
    fname = os.path.join(test.workspace, gwtname + ".obs.csv")
    assert os.path.isfile(fname), f"file not found: {fname}"
    simvals = np.genfromtxt(fname, names=True, delimiter=",", deletechars="")

    sim = test.sims[0]
    gwt = sim.gwt[0]
    mst = gwt.mst
    ist = gwt.ist
    conc = gwt.output.concentration().get_alldata().reshape(200, 300)
    csrb = mst.output.sorbate().get_alldata().reshape(200, 300)
    cim = ist.output.cim().get_alldata().reshape(200, 300)
    cimsrb = ist.output.sorbate().get_alldata().reshape(200, 300)

    # check conc and csrb
    if sorption == "LINEAR":
        csrb_answer = np.where(conc > 0, distcoef * conc, 0)
    if sorption == "FREUNDLICH":
        csrb_answer = np.where(conc > 0, distcoef * conc**sp2, 0)
    if sorption == "LANGMUIR":
        csrb_answer = np.where(
            conc > 0, distcoef * sp2 * conc / (1 + distcoef * conc), 0
        )
    if not np.allclose(csrb[:, 1:], csrb_answer[:, 1:]):
        diff = csrb - csrb_answer
        print("min and max difference")
        print(diff)
        print(diff.min(), diff.max())
        assert False, "csrb not consistent with known answer"

    # check cim and cimsrb
    if sorption == "LINEAR":
        cimsrb_answer = np.where(cim > 0, distcoef * cim, 0)
    if sorption == "FREUNDLICH":
        cimsrb_answer = np.where(cim > 0, distcoef * cim**sp2, 0)
    if sorption == "LANGMUIR":
        cimsrb_answer = np.where(
            cim > 0, distcoef * sp2 * cim / (1 + distcoef * cim), 0
        )
    if not np.allclose(cimsrb[:, 1:], cimsrb_answer[:, 1:]):
        diff = cimsrb - cimsrb_answer
        print("min and max difference")
        print(diff.min(), diff.max())
        assert False, "cimsrb not consistent with known answer"


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

"""
Test basic dispersion for two coupled gwt models.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["dsp01_gwtgwt"]
gdelr = 1.0

# solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def get_gwf_model(sim, gwfname, gwfpath, modelshape):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = gdelr
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xshift,
        yorigin=yshift,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=laytyp,
        k=hk,
        save_specific_discharge=True,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    gwf.set_model_relative_path(gwfpath)
    return gwf


def get_gwt_model(sim, gwtname, gwtpath, modelshape):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0]
    strt = 1.0
    hk = 1.0
    laytyp = 0

    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
    )
    gwt.name_file.save_flows = True

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=gdelr,
        delc=delc,
        top=top,
        botm=botm,
        xorigin=xshift,
        yorigin=yshift,
    )

    # initial conditions
    initial_conc = np.zeros(50)
    if gwtname == "transport1":
        initial_conc[40] = 100.0

    ic = flopy.mf6.ModflowGwtic(gwt, strt=initial_conc)

    # dispersion
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        diffc=100.0,
        alh=0.0,
        alv=0.0,
        ath1=0.0,
        atv=0.0,
        filename=f"{gwtname}.dsp",
    )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    gwt.set_model_relative_path(gwtpath)
    return gwt


def build_models(idx, test):
    # temporal discretization
    nper = 1
    perlen = [5.0]
    nstp = [200]
    tsmult = [1.0]
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=ws, version="mf6", exe_name="mf6", sim_ws=ws)
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=tdis_rc, pname="sim.tdis"
    )

    # solver settings
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    # grid information
    nlay, nrow, ncol = 1, 1, 50

    # Create gwf1 model
    gwf1 = get_gwf_model(sim, "flow1", "flow1", (nlay, nrow, ncol, 0.0, 0.0))

    # Create gwf2 model
    gwf2 = get_gwf_model(sim, "flow2", "flow2", (nlay, nrow, ncol, 50.0 * gdelr, 0.0))

    # gwf-gwf with interface model enabled
    gwfgwf_data = [[(0, 0, ncol - 1), (0, 0, 0), 1, 0.5, 0.5, 1.0, 0.0, 1.0]]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=gwf1.name,
        exgmnameb=gwf2.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="flow1_flow2.gwfgwf",
        dev_interfacemodel_on=True,
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
        filename="flow.ims",
    )
    sim.register_ims_package(imsgwf, [gwf1.name, gwf2.name])

    # Create gwt model
    gwt1 = get_gwt_model(sim, "transport1", "transport1", (nlay, nrow, ncol, 0.0, 0.0))

    # Create gwt model
    gwt2 = get_gwt_model(
        sim, "transport2", "transport2", (nlay, nrow, ncol, 50.0 * gdelr, 0.0)
    )

    # Create GWT GWT exchange
    gwt1gwt2 = flopy.mf6.ModflowGwtgwt(
        sim,
        exgtype="GWT6-GWT6",
        gwfmodelname1=gwf1.name,
        gwfmodelname2=gwf2.name,
        nexg=len(gwfgwf_data),
        exgmnamea=gwt1.name,
        exgmnameb=gwt2.name,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="transport1_transport2.gwtgwt",
    )

    # GWF GWT exchange
    gwfgwt1 = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="flow1",
        exgmnameb="transport1",
        filename="flow1_transport1.gwfgwt",
    )
    gwfgwt2 = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="flow2",
        exgmnameb="transport2",
        filename="flow2_transport2.gwfgwt",
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
        filename="transport.ims",
    )
    sim.register_ims_package(imsgwt, [gwt1.name, gwt2.name])

    return sim, None


def check_output(idx, test):
    gwtname = "transport1"
    fpth = os.path.join(test.workspace, "transport1", f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc1 = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    gwtname = "transport2"
    fpth = os.path.join(test.workspace, "transport2", f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc2 = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # diffusion across both sub-models:
    assert np.all(conc1 > 0.0)
    assert np.all(conc2 > 0.0)

    # no loss of solute
    assert abs(np.sum(conc1) + np.sum(conc2) - 100.0) < 1e-6


@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

"""
Test the advection schemes using multiple flow and transport models
organized from left to right.  Use the well package to inject water
into the first column of cells on the left.  Assign constant head
cells to the last row of the last model.
"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    pytest.param(0, "adv02a_gwtgwt"),
    pytest.param(1, "adv02b_gwtgwt"),
    pytest.param(2, "adv02c_gwtgwt"),
    pytest.param(3, "adv02d_gwtgwt"),
]

scheme = ["upstream", "central", "tvd", "utvd"]
gdelr = 1.0

# solver settings
nouter, ninner = 100, 20
hclose, rclose, relax = 1e-6, 1e-6, 1.0

number_of_models = 5  # this variable is also used to run the parallel test
concentration = 17.0


def get_gwf_model(sim, gwfname, gwfpath, modelshape, chdspd=None, welspd=None):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = gdelr
    delc = 1.0
    delz = 1.0
    top = 1.0
    botm = [top - (k + 1) * delz for k in range(nlay)]
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

    # chd files
    if chdspd is not None:
        chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
            gwf,
            stress_period_data=chdspd,
            save_flows=False,
            pname="CHD-1",
        )

    # wel files
    if welspd is not None:
        wel = flopy.mf6.ModflowGwfwel(
            gwf,
            print_input=True,
            print_flows=True,
            stress_period_data=welspd,
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
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    gwf.set_model_relative_path(gwfpath)
    return gwf


def get_gwt_model(sim, gwtname, gwtpath, modelshape, scheme, sourcerecarray=None):
    nlay, nrow, ncol, xshift, yshift = modelshape
    delr = 1.0
    delc = 1.0
    delz = 1.0
    top = 1.0
    botm = [top - (k + 1) * delz for k in range(nlay)]

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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=concentration)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme)

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    if sourcerecarray is not None:
        ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

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
    perlen = [500.0]
    nstp = [1]
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

    # grid information
    nlay, nrow, ncol = 2, 1, 3

    for imodel in range(number_of_models):
        if imodel == 0:
            welspd = {0: [[(k, 0, 0), 1.0, concentration] for k in range(nlay)]}
        else:
            welspd = None

        if imodel == number_of_models - 1:
            chdspd = {0: [[(k, 0, ncol - 1), 0.0000000] for k in range(nlay)]}
        else:
            chdspd = None

        if imodel == 0:
            sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
        elif imodel == number_of_models - 1:
            sourcerecarray = [()]
        else:
            sourcerecarray = None

        model_name = f"flow{imodel + 1}"
        xshift = imodel * ncol * gdelr
        gwf = get_gwf_model(
            sim,
            model_name,
            model_name,
            (nlay, nrow, ncol, xshift, 0.0),
            chdspd=chdspd,
            welspd=welspd,
        )

        model_name = f"transport{imodel + 1}"
        gwt1 = get_gwt_model(
            sim,
            model_name,
            model_name,
            (nlay, nrow, ncol, xshift, 0.0),
            scheme[idx],
            sourcerecarray=sourcerecarray,
        )

        gwfgwt = flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=f"flow{imodel + 1}",
            exgmnameb=f"transport{imodel + 1}",
            filename=f"flow{imodel + 1}_transport{imodel + 1}.gwfgwt",
        )

    # gwf-gwf with interface model enabled
    for iexg in range(number_of_models - 1):
        gwfgwf_data = [
            [(k, 0, ncol - 1), (0, 0, 0), 1, 0.5, 0.5, 1.0, 0.0, 1.0]
            for k in range(nlay)
        ]
        exgmnamea = f"flow{iexg + 1}"
        exgmnameb = f"flow{iexg + 2}"
        filename = f"{exgmnamea}_{exgmnameb}.gwfgwf"
        gwfgwf = flopy.mf6.ModflowGwfgwf(
            sim,
            exgtype="GWF6-GWF6",
            nexg=len(gwfgwf_data),
            exgmnamea=exgmnamea,
            exgmnameb=exgmnameb,
            exchangedata=gwfgwf_data,
            auxiliary=["ANGLDEGX", "CDIST"],
            filename=filename,
            dev_interfacemodel_on=True,
        )

    # Create GWT GWT exchanges
    for iexg in range(number_of_models - 1):
        exgmnamea = f"transport{iexg + 1}"
        exgmnameb = f"transport{iexg + 2}"
        filename = f"{exgmnamea}_{exgmnameb}.gwtgwt"
        gwtgwt = flopy.mf6.ModflowGwtgwt(
            sim,
            exgtype="GWT6-GWT6",
            gwfmodelname1=f"flow{iexg + 1}",
            gwfmodelname2=f"flow{iexg + 2}",
            adv_scheme=scheme[idx],
            nexg=len(gwfgwf_data),
            exgmnamea=exgmnamea,
            exgmnameb=exgmnameb,
            exchangedata=gwfgwf_data,
            auxiliary=["ANGLDEGX", "CDIST"],
            filename=filename,
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
    sim.register_ims_package(
        imsgwf, [f"flow{nmodel + 1}" for nmodel in range(number_of_models)]
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
    sim.register_ims_package(
        imsgwt,
        [f"transport{nmodel + 1}" for nmodel in range(number_of_models)],
    )

    return sim, None


def check_output(idx, test):
    conclist = []
    for imodel in range(number_of_models):
        gwtname = f"transport{imodel + 1}"
        fpth = pl.Path(test.workspace) / gwtname / f"{gwtname}.ucn"
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
        conclist.append(conc)
    conc_sim = np.hstack(conclist)

    conc_answer = concentration * np.ones(conc_sim.shape)

    diff = conc_sim - conc_answer
    dmax = np.abs(diff).max()

    assert np.allclose(conc_sim, conc_answer), (
        f"Concentrations do not match with known solution. Max diff = {dmax}"
    )


@pytest.mark.parametrize("idx, name", cases)
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

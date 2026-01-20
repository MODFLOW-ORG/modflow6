"""
Test gwt in a subset of the flow grid:

GWT:   XXXXXXXXX     [ transport ]
                           +
GWF: [   flow1   ] + [   flow2   ]

Flow is directed from right to left. The initial concentration equals
zero everywhere except for the rightmost cell. At the end of the run,
all the solute should be in the leftmost cell from the transport model.
(This will be a non-trivial problem in parallel to test the load balance)

"""

import os

import flopy
import pytest
from framework import TestFramework

cases = ["gwt_subset01"]
gdelr = 1.0

# solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def get_gwf_model(sim, gwfname, gwfpath, modelshape, chdspd):
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

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        save_flows=False,
        pname="CHD-1",
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
    delc = 1.0
    top = 1.0
    botm = [0.0]

    strt = ncol * [0.0]
    strt[-1] = 1.0

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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=strt)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="tvd")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=None)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
    )

    gwt.set_model_relative_path(gwtpath)
    return gwt


def build_models(idx, test):
    # temporal discretization
    nper = 1
    perlen = [200.0]
    nstp = [25]
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
    nlay, nrow, ncol = 1, 1, 50

    # Create gwf1 model
    chdspd = {0: [[(0, 0, 0), 0.0]]}
    gwf1 = get_gwf_model(
        sim,
        "flow1",
        "flow1",
        (nlay, nrow, ncol, 0.0, 0.0),
        chdspd=chdspd,
    )

    # Create gwf2 model
    chdspd = {0: [[(0, 0, ncol - 1), 10.0]]}
    gwf2 = get_gwf_model(
        sim,
        "flow2",
        "flow2",
        (nlay, nrow, ncol, 50.0 * gdelr, 0.0),
        chdspd=chdspd,
    )

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
    gwt = get_gwt_model(
        sim,
        "transport",
        "transport",
        (nlay, nrow, ncol, 50.0 * gdelr, 0.0),
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="flow2",
        exgmnameb="transport",
        filename="flow2_transport.gwfgwt",
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
    sim.register_ims_package(imsgwt, [gwt.name])

    return sim, None


def check_output(idx, test):
    gwtname = "transport"

    fpth = os.path.join(test.workspace, gwtname, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # this simply checks if the solute from the right at t = 0
    # has arrived all the way on the left at t = t_end:
    assert conc[0, 0, 0] == pytest.approx(1.0), "concentration should be close to 1.0"


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

"""
Two-dimensional injection of solute into the middle of a square grid.
The test will pass if the results are symmetric.
Based on test_gwt_adv04, this tests the disu package, which
represents a regular MODFLOW grid.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs
from framework import TestFramework

cases = ["disu01a"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 21, 21
    nper = 1
    perlen = [5.0]
    nstp = [200]
    tsmult = [1.0]
    steady = [True]
    delr = np.ones(ncol, dtype=float)
    delc = np.ones(nrow, dtype=float)
    botm = [0.0]
    strt = 1.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 1.0

    top = 1.0
    laytyp = 0

    def get_nn(k, i, j):
        return k * nrow * ncol + i * ncol + j

    # put constant heads all around the box
    chdlist = []
    ib = np.ones((nlay, nrow, ncol), dtype=int)
    ib[:, 1 : nrow - 1, 1 : ncol - 1] = 0
    idloc = np.where(ib > 0)
    for k, i, j in zip(idloc[0], idloc[1], idloc[2]):
        chdlist.append([(get_nn(k, i, j),), 0.0])
    chdspdict = {0: chdlist}

    # injection well with rate and concentration of 1.
    k, i, j = (0, int(nrow / 2), int(ncol / 2))
    w = {0: [[(get_nn(k, i, j),), 1.0, 1.0]]}

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

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

    # use utility to make a disu version of a regular grid
    disu_kwargs = get_disu_kwargs(
        nlay, nrow, ncol, delr, delc, top, botm, return_vertices=True
    )
    disu = flopy.mf6.ModflowGwfdisu(gwf, **disu_kwargs)

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)

    # chd files
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=chdspdict, save_flows=False, pname="CHD-1"
    )

    # wel files
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
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
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

    # use utility to make a disu version of a regular grid
    disu_kwargs = get_disu_kwargs(
        nlay, nrow, ncol, delr, delc, top, botm, return_vertices=True
    )
    disu = flopy.mf6.ModflowGwtdisu(gwt, **disu_kwargs)

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwtname}.adv")

    # dispersion must be off as disu package does not have ANGLDEGX specified
    # dsp = flopy.mf6.ModflowGwtdsp(
    #    gwt,
    #    xt3d_off=True,
    #    diffc=100.0,
    #    alh=0.0,
    #    alv=0.0,
    #    ath1=0.0,
    #    atv=0.0,
    #    filename="{}.dsp".format(gwtname),
    # )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

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
        saverecord=[("CONCENTRATION", "LAST")],
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


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # Check to make sure that the concentrations are symmetric in both the
    # up-down and left-right directions
    conc = conc.reshape((21, 21))
    concud = np.flipud(conc)
    assert np.allclose(concud, conc), (
        "simulated concentrations are not symmetric in up-down direction."
    )

    conclr = np.fliplr(conc)
    assert np.allclose(conclr, conc), (
        "simulated concentrations are not symmetric in left-right direction."
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

"""
Test the SSM package auxiliary variables for specifying source and sink
concentrations for array based recharge.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ssm05"]

nlay, nrow, ncol = 3, 5, 5
idomain_lay0 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 1, 1, 1],
]
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, :, :] = np.array(idomain_lay0)


def build_models(idx, test):
    perlen = [5.0]
    nstp = [5]
    tsmult = [1.0]
    nper = len(perlen)
    delr = 1.0
    delc = 1.0
    top = 4.0
    botm = [3.0, 2.0, 1.0]
    strt = 4.0
    hk = 1.0
    laytyp = 0

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
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
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
        idomain=idomain,
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
    spd = [[(nlay - 1, nrow - 1, ncol - 1), 4.0]]
    chd = flopy.mf6.modflow.ModflowGwfchd(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="CHD-1",
    )

    # array-based rch file
    recharge_rate = np.arange(nrow * ncol).reshape((nrow, ncol))
    recharge_concentration = np.arange(nrow * ncol).reshape((nrow, ncol))
    auxiliary = [("auxvar1", "auxvar2")]
    aux = {0: [[50.0], recharge_concentration]}
    rch2 = flopy.mf6.ModflowGwfrcha(
        gwf,
        print_flows=True,
        recharge=recharge_rate,
        pname="RCH-1",
        filename=f"{gwfname}.rch",
        auxiliary=auxiliary,
        aux=aux,
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

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

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
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt)

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # ssm package
    sourcerecarray = [("RCH-1", "AUX", "AUXVAR2")]
    fileinput = None
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, print_flows=True, sources=sourcerecarray, fileinput=fileinput
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        f"{gwtname}.obs.csv": [
            ("(1-1-1)", "CONCENTRATION", (0, 0, 0)),
            ("(1-5-5)", "CONCENTRATION", (nlay - 1, nrow - 1, ncol - 1)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname=f"{gwtname}.obs",
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


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name

    # load concentration file
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()

    # load transport budget file
    fpth = os.path.join(test.workspace, f"{gwtname}.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    ssmbudall = bobj.get_data(text="SOURCE-SINK MIX")
    times = cobj.get_times()

    print(times)
    for itime, totim in enumerate(times):
        print(f"Checking records for time {totim}")

        # Check records for each of the four recharge packages
        ssmbud = ssmbudall[itime]
        istart = 0
        for irchpak in [0]:
            print(f"  Checking records for recharge package {irchpak + 1}")
            istop = istart + 23

            print(ssmbud[istart:istop])

            print("    Checking id1")
            id1 = ssmbud[istart:istop]["node"]
            id1a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
            assert np.allclose(id1, id1a), f"{id1} /= {id1a}"

            print("    Checking id2")
            id2 = ssmbud[istart:istop]["node2"]
            id2a = id1a
            assert np.allclose(id2, id2a), f"{id2} /= {id2a}"

            print("    Checking q")
            q = ssmbud[istart:istop]["q"]
            qa = [float(a - 1) ** 2 for a in id1a]
            assert np.allclose(q, qa), f"{q} /=\n {qa}"

            istart = istop


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

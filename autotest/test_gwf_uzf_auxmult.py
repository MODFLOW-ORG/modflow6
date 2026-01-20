"""
A simple test of multiple UZF objects within a cell.
The fourth sub-test is designed to fail.  An assert statement
ensure this to be the case.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["uzauxmlt-1to1", "uzauxmlt-2to1", "uzauxmlt-3to1", "uzauxmlt-bad"]
uzarea_correspondences = ("1to1", "2to1", "3to1", "bad")
iuz_cell_dict = {}
cell_iuz_dict = {}

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97

uzarea_data = {
    uzarea_correspondences[0]: {
        "uzf_pkdat": [
            [0, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [1, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [2, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [3, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
        ],
        "auxmultval": 1.0,
    },
    uzarea_correspondences[1]: {
        "uzf_pkdat": [
            [0, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [1, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [2, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [3, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [4, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [5, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [6, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
            [7, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
        ],
        "auxmultval": 0.5,
    },
    uzarea_correspondences[2]: {
        "uzf_pkdat": [
            [0, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [1, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [2, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [3, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [4, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [5, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [6, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [7, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [8, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [9, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
            [10, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
            [11, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
        ],
        "auxmultval": 0.333,
    },
    uzarea_correspondences[3]: {
        "uzf_pkdat": [
            [0, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [1, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [2, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
            [3, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [4, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [5, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
            [6, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [7, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [8, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
            [9, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
            [10, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
            [11, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
        ],
        "auxmultval": 0.334,
    },
}


def build_models(idx, test):
    name = cases[idx]
    uzarea_correspondence = uzarea_correspondences[idx]

    nlay, nrow, ncol = 3, 1, 10
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]

    delr = 100.0
    delc = 100.0
    strt = -9.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON", save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=[-10.0, -20.0, -30.0],
    )

    # initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=1.0e-4,
        k22=1.0e-4,
        k33=1.0e-5,
    )

    # aquifer storage
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=1e-5,
        sy=0.2,
        transient=True,
    )

    # chd files
    chdval = -9.0
    iface = 0
    chdspd = {
        0: [
            [(0, 0, 0), chdval, iface, "object0"],
            [(0, 0, ncol - 1), chdval, iface, "object0"],
        ]
    }
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary="iface",
        boundnames=True,
        print_input=True,
        save_flows=True,
        stress_period_data=chdspd,
    )

    # transient uzf info
    # ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
    uz_pkdat = uzarea_data[uzarea_correspondence]["uzf_pkdat"]

    for itm in uz_pkdat:
        iuz_cell_dict.update({itm[0]: (itm[1][0], itm[1][1], itm[1][2])})
        cell_iuz_dict.update({(itm[1][0], itm[1][1], itm[1][2]): itm[0]})

    finf = 1.0
    extdp = 0.0
    extwc = 0.0
    pet = 0.0
    zero = 0.0
    auxmultval = uzarea_data[uzarea_correspondence]["auxmultval"]
    uzf_spd = {
        0: [
            [i, finf, pet, extdp, extwc, zero, zero, zero, auxmultval]
            for i in np.arange(len(uz_pkdat))
        ]
    }

    flopy.mf6.ModflowGwfuzf(
        gwf,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=True,
        ntrailwaves=7,
        nwavesets=40,
        auxiliary="multiplier",
        auxmultname="multiplier",
        package_convergence_filerecord=f"{name}.UzfConvergence.csv",
        wc_filerecord=f"{name}.wc",
        nuzfcells=len(uz_pkdat),
        packagedata=uz_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{name}.uzf.bud",
        filename=f"{name}.uzf",
    )

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim


def check_output(idx, test):
    print("Running error check")
    name = cases[idx]

    errmsg0 = "flow model should have run successfully but didn't\n"
    errmsg1 = "flow model designed to fail, but seems to have run successfully\n"

    with open(test.workspace / "mfsim.lst", "r") as f:
        lines = f.readlines()
        error_count = 0
        for line in lines:
            if "error report" in line.lower():
                error_count += 1

    if "bad" in name:
        assert error_count > 0, errmsg1
    else:
        assert error_count == 0, errmsg0


@pytest.mark.parametrize("idx,name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    xfail = ["bad" in cases[ct] for ct in np.arange(len(cases))]
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=xfail[idx],
    )
    test.run()

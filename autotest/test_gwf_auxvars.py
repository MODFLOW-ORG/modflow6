import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["aux01"]
auxvar1 = 101.0
auxvar2 = 102.0


def build_models(idx, test):
    nlay, nrow, ncol = 1, 10, 10
    nper = 3
    perlen = [1.0, 1.0, 1.0]
    nstp = [10, 10, 10]
    tsmult = [1.0, 1.0, 1.0]

    lenx = 300.0
    delr = delc = lenx / float(nrow)
    strt = 100.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

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
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
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

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=90.0,
        botm=0.0,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=1.0, k33=0.01)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=0.1,
        steady_state={0: False},
        transient={0: True},
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 100.0])
    chdlist0.append([(0, nrow - 1, ncol - 1), 95.0])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # MAW
    wellbottom = 50.0
    wellrecarray = [[0, 0.1, wellbottom, 100.0, "THIEM", 1, auxvar1, auxvar2]]
    wellconnectionsrecarray = [[0, 0, (0, 5, 5), 100.0, wellbottom, 1.0, 0.1]]
    wellperiodrecarray = [[0, "rate", -0.1]]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord="aux01.maw.bud",
        packagedata=wellrecarray,
        auxiliary=["aux1", "aux2"],
        connectiondata=wellconnectionsrecarray,
        perioddata=wellperiodrecarray,
    )
    # maw.remove()

    # <rno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #       <man> <ncon> <ustrf> <ndv> [<aux(naux)>] [<boundname>]
    packagedata = [
        [
            0,
            (0, 5, ncol - 2),
            delr,
            10.0,
            0.001,
            98.0,
            1.0,
            1.0,
            0.3,
            1,
            1.0,
            0,
            auxvar1,
            auxvar2,
        ],
        [
            1,
            (0, 5, ncol - 1),
            delr,
            10.0,
            0.001,
            97.0,
            1.0,
            1.0,
            0.3,
            1,
            1.0,
            0,
            auxvar1,
            auxvar2,
        ],
    ]
    connectiondata = [[0, -1], [1, 0]]
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_input=True,
        print_stage=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord="aux01.sfr.bud",
        length_conversion=3.281348587,
        time_conversion=86400.0,
        nreaches=len(packagedata),
        packagedata=packagedata,
        auxiliary=["aux1", "aux2"],
        connectiondata=connectiondata,
    )
    # sfr.remove()

    # <ifno> <strt> <nlakeconn> [<aux(naux)>] [<boundname>]
    packagedata = [
        [0, 100.0, 1, auxvar1, auxvar2, "lake1"],
        [1, 100.0, 1, auxvar1, auxvar2, "lake2"],
    ]
    # <ifno> <iconn> <cellid(ncelldim)> <claktype> <bedleak> <belev> <telev> ...
    #        <connlen> <connwidth>
    connectiondata = [
        [0, 0, (0, 1, 1), "vertical", DNODATA, 0.0, 0.0, 0.0, 0.0],
        [1, 0, (0, 2, 2), "vertical", DNODATA, 0.0, 0.0, 0.0, 0.0],
    ]
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        boundnames=True,
        surfdep=1.0,
        print_input=True,
        print_stage=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord="aux01.lak.bud",
        nlakes=len(packagedata),
        packagedata=packagedata,
        auxiliary=["aux1", "aux2"],
        connectiondata=connectiondata,
    )
    # lak.remove()

    # <ifno> <cellid(ncelldim)> <landflag> <ivertcon> <surfdep> <vks> ...
    #        <thtr> <thts> <thti> <eps> [<boundname>]
    packagedata = [
        [0, (0, nrow - 1, 5), 1, -1, 0.1, 0.01, 0.01, 0.1, 0.01, 3.5, "uz1"],
        [1, (0, nrow - 1, 6), 1, -1, 0.1, 0.01, 0.01, 0.1, 0.01, 3.5, "uz1"],
        [2, (0, nrow - 1, 7), 1, -1, 0.1, 0.01, 0.01, 0.1, 0.01, 3.5, "uz1"],
        [3, (0, nrow - 1, 8), 1, -1, 0.1, 0.01, 0.01, 0.1, 0.01, 3.5, "uz1"],
    ]
    # <ifno> <finf> <pet> <extdp> <extwc> <ha> <hroot> <rootact> [<aux(naux)>]
    perioddata = []
    for p in packagedata:
        perioddata.append((p[0], 0.001, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, auxvar1, auxvar2))
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        boundnames=True,
        print_input=True,
        print_flows=True,
        save_flows=True,
        budget_filerecord="aux01.uzf.bud",
        nuzfcells=len(packagedata),
        ntrailwaves=15,
        nwavesets=40,
        packagedata=packagedata,
        auxiliary=["aux1", "aux2"],
        perioddata=perioddata,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    # maw budget aux variables
    fpth = os.path.join(test.workspace, "aux01.maw.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="auxiliary")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)

    # sfr budget aux variables
    fpth = os.path.join(test.workspace, "aux01.sfr.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="auxiliary")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)

    # lak budget aux variables
    fpth = os.path.join(test.workspace, "aux01.maw.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="auxiliary")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)

    # uzf budget aux variables
    fpth = os.path.join(test.workspace, "aux01.uzf.bud")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="auxiliary")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)

    # gwf budget maw aux variables
    fpth = os.path.join(test.workspace, "aux01.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    records = bobj.get_data(text="maw")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)
    records = bobj.get_data(text="sfr")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)
    records = bobj.get_data(text="lak")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)
    records = bobj.get_data(text="uzf")
    for r in records:
        assert np.allclose(r["AUX1"], auxvar1)
        assert np.allclose(r["AUX2"], auxvar2)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

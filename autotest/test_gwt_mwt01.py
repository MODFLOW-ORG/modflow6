"""
Simple 3-layer model with a maw.  Purpose is to test pumping
with concentration being drawn in from edge.  The aquifer
starts with a concentration of zero, but the values grow as the boundary
flows into the aquifer.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mwt_01"]


def build_models(idx, test):
    lx = 5.0
    lz = 3.0
    nlay = 3
    nrow = 1
    ncol = 5
    nper = 1
    delc = 1.0
    delr = lx / ncol
    delz = lz / nlay
    top = [0.0, 0.0, 0.0, 0.0, 0.0]
    botm = list(0 - np.arange(delz, nlay * delz + delz, delz))

    perlen = [0.1]
    nstp = [10]
    kstp = perlen[0] / nstp[0]
    tsmult = [1.0]

    Kh = 20.0
    Kv = 20.0

    steady = [True]
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    single_matrix = False
    nouter, ninner = 700, 300
    hclose, rclose, relax = 1e-8, 1e-6, 0.97

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

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
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
        filename=f"{gwfname}.ims",
    )

    idomain = np.full((nlay, nrow, ncol), 1)
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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=0,
        k=Kh,
        k33=Kv,
    )

    # chd files
    chdlist1 = [
        [(0, 0, 0), 0.0, 100.0],
        [(0, 0, ncol - 1), 0.0, 0.0],
    ]
    chd1 = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdlist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.chd",
    )

    # MAW
    opth = f"{name}.maw.obs"
    wellbottom = -3.0
    wellrecarray = [[0, 0.1, wellbottom, 0.0, "THIEM", 3]]
    wellconnectionsrecarray = [
        [0, 0, (0, 0, 2), 0.0, -1, 1.0, 0.1],
        [0, 1, (1, 0, 2), -1.0, -2, 1.0, 0.1],
        [0, 2, (2, 0, 2), -2.0, -3, 1.0, 0.1],
    ]
    wellperiodrecarray = [[0, "rate", -1.0]]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{gwfname}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        packagedata=wellrecarray,
        connectiondata=wellconnectionsrecarray,
        perioddata=wellperiodrecarray,
        pname="MAW-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
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

    if not single_matrix:
        imsgwt = flopy.mf6.ModflowIms(
            sim,
            print_option="ALL",
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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwtname}.adv")

    # storage
    porosity = 0.30
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")
    # sources
    sourcerecarray = [
        ("CHD-1", "AUX", "CONCENTRATION"),
        # ('WEL-1', 'AUX', 'CONCENTRATION'),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    mwt_obs = {}
    for obstype in [
        "CONCENTRATION",
        "FROM-MVR",
        "STORAGE",
        "CONSTANT",
        "RATE",
        "FW-RATE",
        "RATE-TO-MVR",
        "FW-RATE-TO-MVR",
    ]:
        fname = f"{gwtname}.mwt.obs.{obstype.lower()}.csv"
        ncv = 1
        obs1 = [(f"mwt{i + 1}", obstype, i + 1) for i in range(ncv)]
        obs2 = [(f"bmwt{i + 1}", obstype, f"mymwt{i + 1}") for i in range(ncv)]
        mwt_obs[fname] = obs1 + obs2

    obstype = "MWT"
    fname = f"{gwtname}.mwt.obs.{obstype.lower()}.csv"
    ncv = 1
    nconn = 3
    obs1 = []
    for icv in range(ncv):
        for iconn in range(nconn):
            obs1.append((f"mwt{icv + 1}x{iconn + 1}", obstype, icv + 1, iconn + 1))
    obs2 = [(f"bmwt{i + 1}", obstype, f"mymwt{i + 1}") for i in range(ncv)]
    mwt_obs[fname] = obs1 + obs2

    # append additional obs attributes to obs dictionary
    mwt_obs["digits"] = 15
    mwt_obs["print_input"] = True
    mwt_obs["filename"] = gwtname + ".mwt.obs"

    mwtpackagedata = [
        (0, 0.0, 99.0, 999.0, "mymwt1"),
    ]
    mwtperioddata = [
        (0, "STATUS", "ACTIVE"),
        (0, "CONCENTRATION", 0.0),
    ]
    mwtperioddata = None

    mwt = flopy.mf6.modflow.ModflowGwtmwt(
        gwt,
        boundnames=True,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".mwt.bin",
        budget_filerecord=gwtname + ".mwt.bud",
        packagedata=mwtpackagedata,
        mwtperioddata=mwtperioddata,
        observations=mwt_obs,
        pname="MAW-1",
        auxiliary=["aux1", "aux2"],
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[
            ("CONCENTRATION", "ALL"),
            ("BUDGET", "ALL"),
        ],
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


def check_obs(sim):
    print("checking obs...")
    name = sim.name
    ws = sim.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    # extract mwt concentrations from binary output file
    conc_mwt1 = gwt.mwt.output.concentration().get_alldata().flatten()

    # ensure mwt obs are the same whether specified by
    # boundname or by reach
    csvfiles = gwt.mwt.obs.output.obs_names
    for csvfile in csvfiles:
        if ".mwt.csv" in csvfile:
            continue
        print(f"Checking csv file: {csvfile}")
        conc_ra = gwt.mwt.obs.output.obs(f=csvfile).data
        success = True
        if ".concentration.csv" in csvfile:
            print("Comparing binary concentrations with observed well concentrations.")
            is_same = np.allclose(conc_ra["BMWT1"], conc_mwt1)
            if not is_same:
                success = False
                print(
                    "Binary concentrations do not match with "
                    "observation concentrations for mwt1"
                )
                print(conc_ra["BMWT1"], conc_mwt1)
        # check boundname observations with numeric ID observations
        for icv in range(1):
            # print(f"  Checking reach {imwt + 1}")
            is_same = np.allclose(conc_ra[f"MWT{icv + 1}"], conc_ra[f"BMWT{icv + 1}"])
            if not is_same:
                success = False
                for t, x, y in zip(
                    conc_ra["totim"],
                    conc_ra[f"MWT{icv + 1}"],
                    conc_ra[f"BMWT{icv + 1}"],
                ):
                    print(t, x, y)

    # Sum individual iconn mwt rates and compare with total rate
    csvfile = f"{gwtname}.mwt.obs.mwt.csv"
    print(f"Checking csv file: {csvfile}")
    conc_ra = gwt.mwt.obs.output.obs(f=csvfile).data
    ntimes = conc_ra.shape[0]
    for imwt in range(1):
        connection_sum = np.zeros(ntimes)
        for column_name in conc_ra.dtype.names:
            if f"MWT{icv + 1}X" in column_name:
                connection_sum += conc_ra[column_name]
        is_same = np.allclose(connection_sum, conc_ra[f"BMWT{icv + 1}"])
        if not is_same:
            success = False
            diff = connection_sum - conc_ra[f"BMWT{icv + 1}"]
            print(
                f"Problem with MWT {icv + 1}; "
                f"mindiff {diff.min()} and maxdiff {diff.max()}"
            )

    assert success, "One or more MWT obs checks did not pass"


def check_output(idx, test):
    # ensure mwt concentrations were saved
    name = test.name
    gwtname = "gwt_" + name
    fname = gwtname + ".mwt.bin"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    # ensure gwt concentrations were saved
    fname = gwtname + ".ucn"
    fname = os.path.join(test.workspace, fname)
    assert os.path.isfile(fname)

    check_obs(test)


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

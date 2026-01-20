import os

import flopy
import numpy as np
import pytest
from binary_util import write_budget, write_head
from flopy.utils.gridutil import uniform_flow_field
from framework import TestFramework

cases = ["fmi01a_fc"]
xt3d = [False, True]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 3
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]
    steady = [True]
    delr = 1.0
    delc = 1.0
    top = 1.0
    laytyp = 0
    ss = 0.0
    sy = 0.1
    botm = [0.0]
    strt = 1.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 1.0

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
        idomain=1,
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=10.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt)

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

    # create a heads file with head equal top
    fname = os.path.join(ws, "myhead.hds")
    with open(fname, "wb") as fbin:
        for kstp in range(nstp[0]):
            write_head(fbin, top * np.ones((nrow, ncol)), kstp=kstp + 1)

    # create a budget file
    qx = 0.0
    qy = 0.0
    qz = 0.0
    shape = (nlay, nrow, ncol)
    spdis, flowja = uniform_flow_field(qx, qy, qz, shape)
    dt = np.dtype(
        [
            ("ID1", np.int32),
            ("ID2", np.int32),
            ("FLOW", np.float64),
            ("CONCENTRATION", np.float64),
        ]
    )
    print(flowja)
    # create unbalanced flowja to check flow_imbalance_correction
    # note that residual must be stored in diagonal position
    flowja = np.array([0.01, 0.01, -0.01, 0.0, -0.01, 0.01, 0.01])

    dt = np.dtype(
        [
            ("ID1", np.int32),
            ("ID2", np.int32),
            ("FLOW", np.float64),
            ("SATURATION", np.float64),
        ]
    )
    sat = np.array([(i, i, 0.0, 1.0) for i in range(nlay * nrow * ncol)], dtype=dt)

    fname = os.path.join(ws, "mybudget.bud")
    with open(fname, "wb") as fbin:
        for kstp in range(nstp[0]):
            write_budget(fbin, flowja, kstp=kstp + 1)
            write_budget(fbin, spdis, text="      DATA-SPDIS", imeth=6, kstp=kstp + 1)
            write_budget(fbin, sat, text="        DATA-SAT", imeth=6, kstp=kstp + 1)
    fbin.close()

    # flow model interface
    packagedata = [
        ("GWFBUDGET", "mybudget.bud", None),
        ("GWFHEAD", "myhead.hds", None),
    ]
    fmi = flopy.mf6.ModflowGwtfmi(
        gwt, flow_imbalance_correction=True, packagedata=packagedata
    )

    return sim, None


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()

    # This is the answer to this problem.  Concentration should not change
    cres = [[[10, 10, 10]]]
    cres = np.array(cres)
    errmsg = "simulated concentrations do not match with known solution.\n"
    errmsg += f"cres: {cres}\ncans:{conc}"
    assert np.allclose(cres, conc), errmsg


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

import os

import flopy
import numpy as np
import pytest
from binary_util import write_budget, write_head
from framework import TestFramework

cases = ["dsp01a_fmi", "dsp01b_fmi"]
xt3d = [False, True]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 100
    nper = 1
    perlen = [5.0]
    nstp = [200]
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

    c = {0: [[(0, 0, 0), 0.0000000], [(0, 0, 99), 0.0000000]]}

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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # dispersion
    xt3d_off = not xt3d[idx]
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=xt3d_off,
        diffc=100.0,
        alh=0.0,
        alv=0.0,
        ath1=0.0,
        atv=0.0,
        filename=f"{gwtname}.dsp",
    )

    # constant concentration
    cncs = {0: [[(0, 0, 0), 1.0]]}
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        maxbound=len(cncs),
        stress_period_data=cncs,
        save_flows=False,
        pname="CNC-1",
    )

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # create a heads file with head equal top
    fname = os.path.join(ws, "myhead.hds")
    with open(fname, "wb") as fbin:
        for kstp in range(1):  # nstp[0]):
            write_head(fbin, top * np.ones((nrow, ncol)), kstp=kstp + 1)

    # create a budget file with flows set to zero
    nja = (ncol - 2) * 3 + 2 * 2
    flowja = np.zeros(nja, dtype=np.float64)

    dt = np.dtype(
        [
            ("id1", np.int32),
            ("id2", np.int32),
            ("flow", np.float64),
            ("qx", np.float64),
            ("qy", np.float64),
            ("qz", np.float64),
        ]
    )
    spdis = np.array([(id1, id1, 0.0, 0.0, 0.0, 0.0) for id1 in range(100)], dtype=dt)

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
        for kstp in range(1):  # nstp[0]):
            write_budget(fbin, flowja, kstp=kstp + 1)
            write_budget(fbin, spdis, text="      DATA-SPDIS", imeth=6, kstp=kstp + 1)
            write_budget(fbin, sat, text="        DATA-SAT", imeth=6, kstp=kstp + 1)
    fbin.close()

    # flow model interface
    packagedata = [
        ("GWFBUDGET", "mybudget.bud", None),
        ("GWFHEAD", "myhead.hds", None),
    ]
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=packagedata)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name

    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    try:
        cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
        conc = cobj.get_data()
    except:
        assert False, f'could not load data from "{fpth}"'

    # This is the answer to this problem.  These concentrations are for
    # time step 200.
    cres = [
        [
            [
                1.0,
                0.97472443,
                0.94947431,
                0.92427504,
                0.89915185,
                0.87412972,
                0.84923335,
                0.82448706,
                0.79991471,
                0.77553964,
                0.75138462,
                0.72747174,
                0.70382241,
                0.68045725,
                0.65739608,
                0.63465784,
                0.61226053,
                0.59022124,
                0.56855604,
                0.54727998,
                0.52640705,
                0.50595018,
                0.4859212,
                0.46633085,
                0.44718873,
                0.42850336,
                0.41028211,
                0.39253126,
                0.37525599,
                0.35846038,
                0.34214746,
                0.32631921,
                0.31097658,
                0.29611954,
                0.28174707,
                0.26785727,
                0.2544473,
                0.24151351,
                0.22905142,
                0.21705579,
                0.20552066,
                0.19443937,
                0.18380466,
                0.17360869,
                0.16384304,
                0.15449886,
                0.14556682,
                0.13703721,
                0.12889996,
                0.12114473,
                0.1137609,
                0.10673763,
                0.10006394,
                0.09372869,
                0.08772068,
                0.08202862,
                0.07664126,
                0.07154731,
                0.06673558,
                0.06219493,
                0.05791434,
                0.05388294,
                0.05009,
                0.04652499,
                0.04317758,
                0.04003765,
                0.03709534,
                0.03434103,
                0.03176537,
                0.0293593,
                0.02711402,
                0.02502107,
                0.02307224,
                0.02125967,
                0.01957579,
                0.01801336,
                0.01656542,
                0.01522538,
                0.01398691,
                0.01284404,
                0.01179109,
                0.01082267,
                0.00993375,
                0.00911954,
                0.0083756,
                0.00769775,
                0.00708212,
                0.00652511,
                0.00602341,
                0.00557398,
                0.00517407,
                0.00482116,
                0.00451303,
                0.0042477,
                0.00402344,
                0.00383879,
                0.00369253,
                0.00358368,
                0.00351152,
                0.00347556,
            ]
        ]
    ]
    cres = np.array(cres)
    assert np.allclose(cres, conc), (
        "simulated concentrations do not match with known solution."
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

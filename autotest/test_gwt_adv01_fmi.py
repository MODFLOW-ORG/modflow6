"""
Test the advection schemes in the gwt advection package for a one-dimensional
model grid of square cells.
"""

import os

import flopy
import numpy as np
import pytest
from binary_util import write_budget, write_head
from flopy.utils.gridutil import uniform_flow_field
from framework import TestFramework

cases = ["adv01a_fmi", "adv01b_fmi", "adv01c_fmi"]
scheme = ["upstream", "central", "tvd"]


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
    botm = [0.0]
    strt = 1.0
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

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme[idx], filename=f"{gwtname}.adv")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # create a heads file with head equal top
    fname = os.path.join(ws, "myhead.hds")
    with open(fname, "wb") as fbin:
        for kstp in range(1):  # nstp[0]):
            write_head(fbin, top * np.ones((nrow, ncol)), kstp=kstp + 1)

    # create a budget file
    qx = 1.0
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
    wel = np.array([(0 + 1, 0 + 1, 1.0, 1.0)], dtype=dt)
    chd = np.array([(ncol - 1 + 1, 0 + 1, -1.0, 0.0)], dtype=dt)

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
            write_budget(
                fbin,
                wel,
                text="             WEL",
                imeth=6,
                text2id2="           WEL-1",
                kstp=kstp + 1,
            )
            write_budget(
                fbin,
                chd,
                text="             CHD",
                imeth=6,
                text2id2="           CHD-1",
                kstp=kstp + 1,
            )
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

    obs_data = {
        "conc_obs.csv": [
            ("(1-1-10)", "CONCENTRATION", (0, 0, 9)),
            ("(1-1-50)", "CONCENTRATION", (0, 0, 49)),
        ],
        "flow_obs.csv": [
            ("c10-c11", "FLOW-JA-FACE", (0, 0, 9), (0, 0, 10)),
            ("c50-c51", "FLOW-JA-FACE", (0, 0, 49), (0, 0, 50)),
            ("c99-c100", "FLOW-JA-FACE", (0, 0, 98), (0, 0, 99)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname="conc_obs",
        filename=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
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
    cres1 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999997e-01,
                9.99999991e-01,
                9.99999971e-01,
                9.99999914e-01,
                9.99999761e-01,
                9.99999372e-01,
                9.99998435e-01,
                9.99996286e-01,
                9.99991577e-01,
                9.99981712e-01,
                9.99961893e-01,
                9.99923632e-01,
                9.99852532e-01,
                9.99725120e-01,
                9.99504599e-01,
                9.99135431e-01,
                9.98536850e-01,
                9.97595635e-01,
                9.96158712e-01,
                9.94026505e-01,
                9.90948130e-01,
                9.86619748e-01,
                9.80687319e-01,
                9.72754814e-01,
                9.62398489e-01,
                9.49187176e-01,
                9.32707801e-01,
                9.12594513e-01,
                8.88559134e-01,
                8.60420154e-01,
                8.28127324e-01,
                7.91779115e-01,
                7.51630867e-01,
                7.08092322e-01,
                6.61714306e-01,
                6.13165405e-01,
                5.63200494e-01,
                5.12623768e-01,
                4.62249349e-01,
                4.12862664e-01,
                3.65185517e-01,
                3.19847250e-01,
                2.77363614e-01,
                2.38124183e-01,
                2.02388273e-01,
                1.70288648e-01,
                1.41841739e-01,
                1.16962748e-01,
                9.54838854e-02,
                7.71740354e-02,
                6.17583229e-02,
                4.89363652e-02,
                3.83983188e-02,
                2.98381826e-02,
                2.29641338e-02,
                1.75059339e-02,
                1.32196416e-02,
                9.89000005e-03,
                7.33093269e-03,
                5.38459977e-03,
                3.91944360e-03,
                2.82760119e-03,
                2.02199855e-03,
                1.43337156e-03,
                1.00739149e-03,
                7.02013580e-04,
                4.85116958e-04,
                3.32465664e-04,
                2.25991387e-04,
                1.52379541e-04,
                1.01928496e-04,
                6.76460984e-05,
                4.45462926e-05,
                2.91101871e-05,
                1.88792800e-05,
                1.21527525e-05,
                7.76522212e-06,
                4.92565188e-06,
                3.10201677e-06,
                1.93969988e-06,
                1.20440812e-06,
                7.42676511e-07,
                4.54831064e-07,
                2.76669882e-07,
                1.67174989e-07,
                1.00349240e-07,
                5.98446532e-08,
                3.54600737e-08,
            ]
        ]
    ]
    cres1 = np.array(cres1)

    cres2 = [
        [
            [
                9.99996617e-01,
                1.00001184e00,
                1.00000294e00,
                9.99972914e-01,
                9.99992627e-01,
                1.00004237e00,
                1.00002081e00,
                9.99945149e-01,
                9.99952654e-01,
                1.00005669e00,
                1.00008810e00,
                9.99966402e-01,
                9.99865541e-01,
                9.99967791e-01,
                1.00015792e00,
                1.00014755e00,
                9.99895530e-01,
                9.99724106e-01,
                9.99916592e-01,
                1.00029941e00,
                1.00038455e00,
                9.99960678e-01,
                9.99433053e-01,
                9.99453350e-01,
                1.00018163e00,
                1.00097923e00,
                1.00093550e00,
                9.99790199e-01,
                9.98371554e-01,
                9.98054584e-01,
                9.99598363e-01,
                1.00229288e00,
                1.00416575e00,
                1.00323035e00,
                9.98995210e-01,
                9.93234271e-01,
                9.89448228e-01,
                9.91206357e-01,
                1.00016889e00,
                1.01473298e00,
                1.02990960e00,
                1.03846239e00,
                1.03282855e00,
                1.00710727e00,
                9.58480908e-01,
                8.87726436e-01,
                7.98820097e-01,
                6.97900399e-01,
                5.91969549e-01,
                4.87686471e-01,
                3.90487541e-01,
                3.04127133e-01,
                2.30608327e-01,
                1.70400015e-01,
                1.22812141e-01,
                8.64138068e-02,
                5.94120233e-02,
                3.99463958e-02,
                2.62868102e-02,
                1.69426845e-02,
                1.07033555e-02,
                6.63198283e-03,
                4.03300421e-03,
                2.40844447e-03,
                1.41323306e-03,
                8.15254552e-04,
                4.62589305e-04,
                2.58303233e-04,
                1.42001900e-04,
                7.68911977e-05,
                4.10256980e-05,
                2.15775541e-05,
                1.11912143e-05,
                5.72578796e-06,
                2.89083689e-06,
                1.44073067e-06,
                7.09001789e-07,
                3.44624235e-07,
                1.65501321e-07,
                7.85475047e-08,
                3.68512253e-08,
                1.70949923e-08,
                7.84310280e-09,
                3.55966819e-09,
                1.59856594e-09,
                7.10467596e-10,
                3.12565151e-10,
                1.36146377e-10,
                5.87252052e-11,
                2.50886169e-11,
                1.06179506e-11,
                4.45237718e-12,
                1.85013624e-12,
                7.61982955e-13,
                3.11095972e-13,
                1.25908830e-13,
                5.05704707e-14,
                2.00370648e-14,
                8.15003576e-15,
                2.57563506e-15,
            ]
        ]
    ]
    cres2 = np.array(cres2)

    cres3 = [
        [
            [
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                1.00000000e00,
                9.99999999e-01,
                9.99999997e-01,
                9.99999991e-01,
                9.99999975e-01,
                9.99999926e-01,
                9.99999789e-01,
                9.99999407e-01,
                9.99998374e-01,
                9.99995665e-01,
                9.99988785e-01,
                9.99971918e-01,
                9.99932078e-01,
                9.99841550e-01,
                9.99643930e-01,
                9.99229970e-01,
                9.98398720e-01,
                9.96800070e-01,
                9.93857995e-01,
                9.88681096e-01,
                9.79978744e-01,
                9.66015902e-01,
                9.44652308e-01,
                9.13514114e-01,
                8.70328697e-01,
                8.13410724e-01,
                7.42224214e-01,
                6.57879960e-01,
                5.63390876e-01,
                4.63530320e-01,
                3.64233335e-01,
                2.71628522e-01,
                1.90935412e-01,
                1.25541007e-01,
                7.65316248e-02,
                4.28052252e-02,
                2.16851758e-02,
                9.78976172e-03,
                3.85613094e-03,
                1.28872611e-03,
                3.52070089e-04,
                7.49188445e-05,
                1.17688715e-05,
                1.33952025e-06,
                1.08174095e-07,
                -4.82019087e-08,
                -5.67180537e-08,
                -4.65251289e-08,
                -3.25511455e-08,
                -1.94644548e-08,
                -9.78876693e-09,
                -4.07380361e-09,
                -1.38097809e-09,
                -3.72934181e-10,
                -7.83508455e-11,
                -1.26040926e-11,
                -1.48260453e-12,
                4.10392230e-14,
                2.44993743e-13,
                2.46295025e-13,
                1.90964563e-13,
                1.03476379e-13,
                3.96502895e-14,
                1.04500247e-14,
                2.00830327e-15,
                4.70831032e-16,
                3.38440506e-16,
                2.49848438e-16,
                1.83245111e-16,
                1.32361223e-16,
                9.39406563e-17,
                6.54891851e-17,
                4.48667613e-17,
                3.02333440e-17,
                2.00567815e-17,
                1.31110206e-17,
                8.45177289e-18,
                5.37610069e-18,
                3.37597383e-18,
            ]
        ]
    ]
    cres3 = np.array(cres3)

    creslist = [cres1, cres2, cres3]

    assert np.allclose(creslist[idx], conc), (
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

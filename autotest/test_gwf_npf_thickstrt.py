import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "gwf_npf_thickstrt01",  # icelltype=0
    "gwf_npf_thickstrt02",  # icelltype=0, using thickstrt, but it has no effect
    "gwf_npf_thickstrt03",  # icelltype=-1, using thickstrt and strt = 5.
    "gwf_npf_thickstrt04",  # icelltype=1, no thickstrt and strt = 5.
    "gwf_npf_thickstrt05",  # icelltype=-1, no thickstrt and strt = 5.
    "gwf_npf_thickstrt06",  # icelltype=0, no thickstrt, has hfb
    "gwf_npf_thickstrt07",  # icelltype=-1, using thickstrt, has hfb
    "gwf_npf_thickstrt08",  # icelltype=1, no thickstrt, has hfb
    "gwf_npf_thickstrt09",  # icelltype=-1, no thickstrt, has hfb
]
thickstrt = [False, True, True, False, False, False, True, False, False]
icelltype = [0, 0, -1, 1, -1, 0, -1, 1, -1]
hfb_on = [False, False, False, False, False, True, True, True, True]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 6
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = 10.0
    botm = [0.0]
    strt = 5.0
    hk = 1.0

    c = {0: [[(0, 0, 0), 6.0], [(0, 0, ncol - 1), 4.0]]}

    nouter, ninner = 10, 5
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "flow"

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )
    gwf.name_file.save_flows = True

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
        filename=f"{name}.ims",
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
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # npf
    thickstrt_option = thickstrt[idx]
    ict = icelltype[idx]
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, thickstrt=thickstrt_option, icelltype=ict, k=hk, k33=hk
    )

    if hfb_on[idx]:
        hfb = flopy.mf6.ModflowGwfhfb(
            gwf,
            print_input=True,
            maxhfb=1,
            stress_period_data=[((0, 0, 2), (0, 0, 3), 1.0e-4)],
        )

    # chd files
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        print_flows=True,
        pname="CHD-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    name = "flow"

    fpth = os.path.join(test.workspace, f"{name}.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    head = hobj.get_data().flatten()

    # This is the answer to this problem.
    answer_linear = np.linspace(6, 4, 6)

    answer_water_table = (6.0, 5.65716, 5.29206, 4.89969, 4.47276, 4.0)
    answer_water_table = np.array(answer_water_table)

    answer_confined_hfb = (6.0, 5.9998, 5.9996, 4.0004, 4.0002, 4.0)
    answer_confined_hfb = np.array(answer_confined_hfb)

    answer_confined_thickstart_hfb = (
        6.0,
        5.9996004,
        5.9992008,
        4.0007992,
        4.0003996,
        4.0,
    )
    answer_confined_thickstart_hfb = np.array(answer_confined_thickstart_hfb)

    answer_unconfined_hfb = (6.0, 5.99983342, 5.99966683, 4.00049971, 4.00024986, 4.0)
    answer_unconfined_hfb = np.array(answer_unconfined_hfb)

    answer_dict = {
        0: answer_linear,
        1: answer_linear,
        2: answer_linear,
        3: answer_water_table,
        4: answer_water_table,
        5: answer_confined_hfb,
        6: answer_confined_thickstart_hfb,
        7: answer_unconfined_hfb,
        8: answer_unconfined_hfb,
    }

    hres = answer_dict[idx]
    assert np.allclose(hres, head), "simulated head do not match with known solution."

    fpth = os.path.join(test.workspace, f"{name}.cbc")
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")
    q_simulated_inflow = cobj.get_data(idx=1)[0]["q"][0]
    q_answer_dict = {
        0: 4.0,
        1: 4.0,
        2: 2.0,
        3: 1.9965396769631871,
        4: 1.9965396769631871,
        5: 1.9990e-03,
        6: 1.9980e-03,
        7: 9.9949e-04,
        8: 9.9949e-04,
    }
    q_answer = q_answer_dict[idx]
    assert np.allclose(q_answer, q_simulated_inflow), (
        "simulated flow does not match with known solution."
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

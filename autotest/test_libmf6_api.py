"""
Test the MODFLOW API with a basic 1x1x5 model:

  [ BC | 1 | 2 | 3 | BC ]

"""

import os

import flopy
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libmf6_api"]


def get_model(dir):
    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose = 10e-9, 1e-3

    # model spatial discretization
    nlay = 1
    nrow = 1
    ncol = 5

    # cell spacing
    delr = 1.0
    delc = 1.0

    # top/bot of the aquifer
    tops = [0.0, -1.0]

    # hydraulic conductivity
    k11 = 1.0

    # boundary stress period data
    h_left = 0.0
    h_right = 5.0

    # initial head
    h_start = 0.0

    sim = flopy.mf6.MFSimulation(
        sim_name="sim", version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
    )

    # submodel on the left:
    chd_data = [[(0, 0, 0), h_left], [(0, 0, ncol - 1), h_right]]
    chd_spd = {0: chd_data}

    gwf = flopy.mf6.ModflowGwf(sim, modelname="model", save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="model.hds",
        budget_filerecord="model.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws)

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    sim_compare = get_model(ws)

    return sim, sim_compare


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    if model_ws is None:
        model_ws = "."
    output_file_path = os.path.join(model_ws, "mfsim.stdout")

    try:
        mf6 = ModflowApi(exe, working_directory=model_ws)
    except Exception as e:
        print("Failed to load " + str(exe))
        print("with message: " + str(e))
        return False, open(output_file_path).readlines()

    # initialize the model
    try:
        mf6.initialize()
    except:
        return False, open(output_file_path).readlines()

    # testing API
    comp_str = mf6.get_component_name()
    version_str = mf6.get_version()
    print(f"Loaded {comp_str} with version {version_str}")

    # time loop
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        try:
            mf6.update()
        except:
            return False, open(output_file_path).readlines()
        current_time = mf6.get_current_time()

    # finish
    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    # cleanup and return
    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()

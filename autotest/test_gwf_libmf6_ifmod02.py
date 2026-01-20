"""
Test the interface model approach for an inhomogeneous coupling
of three gwf models using the API. One exchange will have XT3D
enabled (Exg1) and the other one (Exg2) doesn't. And the top-left
model's NPF will be configured with XT3D as well. The numbers
in the figure are the global (solution level) indices:

      'top-left'       'top-right'
                   |
    01 02 03 04 05 |  21 22 23 24 25
                  Exg1
    06 07 08 09 10 |  26 27 28 29 30
                   |
    -----Exg2----

    11 12 13 14 15
                          y
    16 17 18 19 20        |
                          |
      'bottom'            0 ---- x

Note that this global numbering is determined by the order in which
the models are added to the simulation below.

From this configuration, we will assert properties such as:

- model 'top-left' does not have one but two interface models
- the grid sizes for it are
    o 12: {3,4,5,21,22,8,9,10,26,27,14,15} for the interface model based on Exg1
    o 10: {6,7,8,9,10,11,12,13,14,15}      for the interface model based on Exg2
  and the grid for the interface model for Exg1 thus contains cells
  from model 'bottom'
- the coefficients for the internal connection between cell 9 and cell 10
  are calculated by the model itself (as it does not have XT3D enabled)
- ...

"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_ifmod02"]

# global convenience...
name_tl = "topleft"
name_tr = "topright"
name_bl = "bottomleft"


def get_model(dir, name):
    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = 10e-9, 1e-3, 0.97

    # model spatial discretization
    nlay = 1
    ncol = 5
    nrow = 2

    # cell spacing
    delr = 8.0
    delc = 13.0
    area = delr * delc

    # top/bot of the aquifer
    tops = [0.0, -5.0]

    # hydraulic conductivity
    k11 = 10.0

    # boundary stress period data
    h_left = -1.0
    h_right = -2.0

    # initial head
    h_start = -2.1

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=dir,
        memory_print_option="all",
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

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
        relaxation_factor=relax,
    )

    # boundary for submodels on the left:
    left_chd = [
        [(ilay, irow, 0), h_left] for irow in range(nrow) for ilay in range(nlay)
    ]
    chd_spd_left = {0: left_chd}

    # boundary for submodel on the right:
    right_chd = [
        [(ilay, irow, ncol - 1), h_right]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    chd_spd_right = {0: right_chd}

    # --------------------------------------
    # top-left model
    # --------------------------------------
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name_tl, save_flows=True)
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
        xt3doptions=True,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_left)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name_tl}.hds",
        budget_filerecord=f"{name_tl}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # --------------------------------------
    # bottom-left model
    # --------------------------------------
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name_bl, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        xorigin=0.0,
        yorigin=-2 * delc,
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
        xt3doptions=False,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_left)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name_bl}.hds",
        budget_filerecord=f"{name_bl}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # --------------------------------------
    # top-right model
    # --------------------------------------
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name_tr, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        xorigin=5 * delr,
        yorigin=0.0,
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
        xt3doptions=False,
        save_flows=True,
        icelltype=0,
        k=k11,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_right)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name_tr}.hds",
        budget_filerecord=f"{name_tr}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # exchange between top-left and top-right (Exg1)
    angldegx = 0.0
    cdist = delr
    gwfgwf_data_exg1 = [
        [
            (ilay, irow, ncol - 1),
            (ilay, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    exg1 = flopy.mf6.ModflowGwfgwf(
        sim,
        filename=name_tl + "-" + name_tr + ".gwfgwf",
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data_exg1),
        exgmnamea=name_tl,
        exgmnameb=name_tr,
        exchangedata=gwfgwf_data_exg1,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=True,  # this will activate the interface model
    )

    # exchange between top-left and bottom-left (Exg2)
    angldegx = 270.0
    cdist = delc
    gwfgwf_data_exg2 = [
        [
            (ilay, nrow - 1, icol),
            (ilay, 0, icol),
            1,
            delc / 2.0,
            delc / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for icol in range(ncol)
        for ilay in range(nlay)
    ]
    exg2 = flopy.mf6.ModflowGwfgwf(
        sim,
        filename=name_tl + "-" + name_bl + ".gwfgwf",
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data_exg2),
        exgmnamea=name_tl,
        exgmnameb=name_bl,
        exchangedata=gwfgwf_data_exg2,
        auxiliary=["ANGLDEGX", "CDIST"],
        dev_interfacemodel_on=True,  # by default no interface model, we force it here
    )

    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name)

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    sim_compare = get_model(ws, name)

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

    # test the interface models
    check_interface_models(mf6)

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


def check_interface_models(mf6):
    # interface model for model top-left from exchange=Exg1
    exg1_id = 1
    ifm_topleft_1 = f"GWFIM1_{exg1_id}"
    gfc_topleft_1 = f"GWFCON1_{exg1_id}"
    exg2_id = 2
    ifm_topleft_2 = f"GWFIM1_{exg2_id}"
    gfc_topleft_2 = f"GWFCON1_{exg2_id}"

    # interface model nr.1
    addr = mf6.get_var_address("IDXTOGLOBALIDX", gfc_topleft_1, "GC")
    idxToGlobalIdx_1 = mf6.get_value_ptr(addr)
    assert np.size(idxToGlobalIdx_1) == 12
    assert np.array_equal(idxToGlobalIdx_1, [3, 4, 5, 21, 22, 8, 9, 10, 26, 27, 14, 15])
    addr = mf6.get_var_address("ia", gfc_topleft_1, "gc/con")
    ia_1 = mf6.get_value_ptr(addr)
    addr = mf6.get_var_address("ja", gfc_topleft_1, "gc/con")
    ja_1 = mf6.get_value_ptr(addr)
    addr = mf6.get_var_address("mask", gfc_topleft_1, "gc/con")
    mask_1 = mf6.get_value_ptr(addr)

    # we have the following connections to calculate,
    # check if the mask is set correctly
    connections = [
        (2, 3),
        (3, 2),
        (3, 4),
        (4, 3),
        (3, 8),
        (8, 3),
        (7, 8),
        (8, 7),
        (8, 9),
        (9, 8),
    ]
    for n in range(1, 13):
        for ipos in range(ia_1[n - 1] + 1, ia_1[n]):
            m = ja_1[ipos - 1]
            mask_nm = mask_1[ipos - 1]
            if mask_nm > 0:
                assert (n, m) in connections
                connections.remove((n, m))

    # so now the list should be empty again
    assert len(connections) == 0

    # interface model nr.2
    addr = mf6.get_var_address("IDXTOGLOBALIDX", gfc_topleft_2, "GC")
    idxToGlobalIdx_2 = mf6.get_value_ptr(addr)
    assert np.size(idxToGlobalIdx_2) == 10
    assert np.array_equal(idxToGlobalIdx_2, [6, 7, 8, 9, 10, 11, 12, 13, 14, 15])


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()

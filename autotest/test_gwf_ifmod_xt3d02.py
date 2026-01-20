"""
Test the interface model approach.
It compares the result of a single, strongly anisotropic model
with XT3D enabled to the equivalent case where the domain is
decomposed and joined by a GWF-GWF exchange with XT3D applied.

       'refmodel'              'leftmodel'     'rightmodel'

   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1    VS    1 1 1 1 1   +   1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1
   1 1 1 1 1 1 1 1 1 1          1 1 1 1 1       1 1 1 1 1

The head values should always be identical. All models are
part of the same solution for convenience.
In addition, a check on the x,y,z components of specific discharge
is present. The values of the left submodel are compared to
the left part of the full model, and similar for right: they
should be identical. Finally, the budget error is checked.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ifmod_xt3d02"]
mname_ref = "refmodel"
mname_left = "leftmodel"
mname_right = "rightmodel"
hclose_check = 1e-9
max_inner_it = 300
useXT3D = True


def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, max_inner_it
    hclose, rclose, relax = hclose_check, 1e-3, 0.97

    # model spatial discretization
    nlay = 1
    ncol = 10
    ncol_left = 5
    ncol_right = 5
    nrow = 10

    # cell spacing
    delr = 10.0
    delc = 10.0
    area = delr * delc

    # shift (hor. and vert.)
    shift_some_x = -20 * delr  # avoids overlap
    shift_x = 5 * delr
    shift_y = 0.0

    # top/bot of the aquifer
    tops = [0.0, -5.0]

    # hydraulic conductivity
    k11 = 10.0
    k22 = 0.1
    k_angle = 45.0

    # boundary stress period data
    h_left = -2.0
    h_right = -2.0

    # initial head
    h_start = -2.0

    # well
    well_id = (0, 4, 4)
    well_rate = -1.0

    # This creates the single model, for reference:
    left_chd = [[(0, irow, 0), h_left] for irow in range(nrow)]
    right_chd = [[(0, irow, ncol - 1), h_right] for irow in range(nrow)]
    chd_data = left_chd + right_chd
    chd_spd = {0: chd_data}

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
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

    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_ref, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        xorigin=shift_some_x,
        yorigin=0.0,
        top=tops[0],
        botm=tops[1:],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=useXT3D,
        save_flows=True,
        icelltype=0,
        k=k11,
        k22=k22,
        angle1=k_angle,
    )

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # well
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=[[well_id, well_rate]],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_ref}.hds",
        budget_filerecord=f"{mname_ref}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # Now create two coupled models with the interface model enabled,
    # to be stored in the same solution as the reference model

    # submodel on the left:
    left_chd = [[(0, irow, 0), h_left] for irow in range(nrow)]
    chd_spd_left = {0: left_chd}

    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_left, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_left,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=useXT3D,
        save_flows=True,
        icelltype=0,
        k=k11,
        k22=k22,
        angle1=k_angle,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_left)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_left}.hds",
        budget_filerecord=f"{mname_left}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=[[well_id, well_rate]],
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
    )

    # submodel on the right:
    right_chd = [[(0, irow, ncol_right - 1), h_right] for irow in range(nrow)]
    chd_spd_right = {0: right_chd}

    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_right, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol_right,
        delr=delr,
        delc=delc,
        xorigin=shift_x,
        yorigin=shift_y,
        top=tops[0],
        botm=tops[1:],
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=useXT3D,
        save_flows=True,
        icelltype=0,
        k=k11,
        k22=k22,
        angle1=k_angle,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_right)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_right}.hds",
        budget_filerecord=f"{mname_right}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # exchangedata
    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (0, irow, ncol_left - 1),
            (0, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for irow in range(nrow)
    ]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=mname_left,
        exgmnameb=mname_right,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        xt3d=useXT3D,
    )

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def qxqyqz(fname, nlay, nrow, ncol):
    nodes = nlay * nrow * ncol
    cbb = flopy.utils.CellBudgetFile(fname, precision="double")
    spdis = cbb.get_data(text="DATA-SPDIS")[0]
    qx = np.ones((nodes), dtype=float) * 1.0e30
    qy = np.ones((nodes), dtype=float) * 1.0e30
    qz = np.ones((nodes), dtype=float) * 1.0e30
    n0 = spdis["node"] - 1
    qx[n0] = spdis["qx"]
    qy[n0] = spdis["qy"]
    qz[n0] = spdis["qz"]
    qx = qx.reshape(nlay, nrow, ncol)
    qy = qy.reshape(nlay, nrow, ncol)
    qz = qz.reshape(nlay, nrow, ncol)
    qx = np.ma.masked_equal(qx, 1.0e30)
    qy = np.ma.masked_equal(qy, 1.0e30)
    qz = np.ma.masked_equal(qz, 1.0e30)
    return qx, qy, qz


def check_output(idx, test):
    print("comparing heads and spec. discharge to single model reference...")

    fpth = os.path.join(test.workspace, f"{mname_ref}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_ref}.cbc")
    nlay, nrow, ncol = heads.shape
    qxb, qyb, qzb = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{mname_left}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_left = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_left}.cbc")
    nlay, nrow, ncol = heads_left.shape
    qxb_left, qyb_left, qzb_left = qxqyqz(fpth, nlay, nrow, ncol)

    fpth = os.path.join(test.workspace, f"{mname_right}.hds")
    hds = flopy.utils.HeadFile(fpth)
    heads_right = hds.get_data()
    fpth = os.path.join(test.workspace, f"{mname_right}.cbc")
    nlay, nrow, ncol = heads_right.shape
    qxb_right, qyb_right, qzb_right = qxqyqz(fpth, nlay, nrow, ncol)

    heads_2models = np.append(heads_left[0], heads_right[0], axis=1)

    # compare heads
    maxdiff = np.amax(abs(heads - heads_2models))
    assert maxdiff < 10 * hclose_check, (
        f"Max. head diff. {maxdiff} should \
                     be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_x left
    maxdiff = np.amax(abs(qxb[:, :, 0:5] - qxb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (x) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_y left
    maxdiff = np.amax(abs(qyb[:, :, 0:5] - qyb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (y) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_z left
    maxdiff = np.amax(abs(qzb[:, :, 0:5] - qzb_left))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (z) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_x right
    maxdiff = np.amax(abs(qxb[:, :, 5:] - qxb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (x) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_y right
    maxdiff = np.amax(abs(qyb[:, :, 5:] - qyb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (y) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # compare spdis_z right
    maxdiff = np.amax(abs(qzb[:, :, 5:] - qzb_right))
    assert maxdiff < 10 * hclose_check, (
        f"Max. diff. in spec. discharge (z) {maxdiff} \
                     should be within solver tolerance (x10): {10 * hclose_check}"
    )

    # check budget error from .lst file
    for mname in [mname_ref, mname_left, mname_right]:
        fpth = os.path.join(test.workspace, f"{mname}.lst")
        for line in open(fpth):
            if line.lstrip().startswith("PERCENT"):
                cumul_balance_error = float(line.split()[3])
                assert abs(cumul_balance_error) < 0.00001, (
                    f"Cumulative balance error = {cumul_balance_error} for {mname}, "
                    "should equal 0.0"
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

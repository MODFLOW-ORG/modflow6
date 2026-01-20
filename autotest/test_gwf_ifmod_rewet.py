"""
General test for the interface model approach.
It compares the result of a single reference model
to the equivalent case where the domain is decomposed
and joined by a GWF-GWF exchange.

In this case we test rewetting, which is also enabled in
the interface model and should give identical results.

period 1: The first stress period we start almost dry and have the
          model fill up.
period 2: The BC on the left is lowered such that a part of the top
          layer dries. To test the interface, the value is chosen such
          that the boundary cell on the left is DRY and the one on the
          right isn't.

                 'refmodel'               'leftmodel'    'rightmodel'

   layer 1:  1 . . . . . . . 1            1 . . . . 1       1 . . 1
   layer 2:  1 . . . . . . . 1     VS     1 . . . . 1   +   1 . . 1
   layer 3:  1 . . . . . . . 1            1 . . . . 1       1 . . 1

We assert equality on the head values. All models are part of the same
solution for convenience. Finally, the budget error is checked.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ifmod_rewet01"]

# some global convenience...:
# model names
mname_ref = "refmodel"
mname_left = "leftmodel"
mname_right = "rightmodel"

# solver criterion
hclose_check = 1e-9
max_inner_it = 300
nper = 2

# model spatial discretization
nlay = 3
ncol = 15
ncol_left = 10
ncol_right = 5
nrow = 1

lenx = 15.0 * 500.0
leny = 10.0 * 500.0
delr = lenx / float(ncol)
delc = leny / float(nrow)
area = delr * delc

# shift (hor. and vert.)
shift_some_x = -25 * delr  # avoids overlap
shift_x = ncol_left * delr
shift_y = 0.0

# top/bot of the aquifer
tops = [150.0, 50.0, 0.0, -50.0]

# hydraulic conductivity
hk = 10.0

# boundary stress period data for period 1 and 2
h_left = [150.0, 20.0]
h_right = 60.0

# initial head
h_start = -40.0


# head boundaries
lchd1 = [
    [(ilay, irow, 0), h_left[0]]
    for ilay in range(nlay)
    for irow in range(nrow)
    if h_left[0] > tops[ilay + 1]
]
rchd = [
    [(ilay, irow, ncol - 1), h_right]
    for ilay in range(nlay)
    for irow in range(nrow)
    if h_right > tops[ilay + 1]
]
rchd_right = [
    [(ilay, irow, ncol_right - 1), h_right]
    for ilay in range(nlay)
    for irow in range(nrow)
    if h_right > tops[ilay + 1]
]
chd1 = lchd1 + rchd

chd_spd = {0: chd1}
chd_spd_left = {0: lchd1}
chd_spd_right = {0: rchd_right}

lchd2 = [
    [(ilay, irow, 0), h_left[1]]
    for ilay in range(nlay)
    for irow in range(nrow)
    if h_left[1] > tops[ilay + 1]
]
chd_spd[1] = lchd2 + rchd
chd_spd_left[1] = lchd2
chd_spd_right[1] = rchd_right

# rewetting
rewet_record = [("WETFCT", 1.0, "IWETIT", 1, "IHDWET", 1)]
wetdry = -0.001


def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, max_inner_it
    hclose, rclose, relax = hclose_check, 1e-3, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
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
        filename="gwf.ims",
    )

    # the full gwf model as a reference
    add_refmodel(sim)

    # now add two coupled models with the interface model enabled,
    # to be stored in the same solution as the reference model
    add_leftmodel(sim)
    add_rightmodel(sim)
    add_gwfexchange(sim)

    return sim


def add_refmodel(sim):
    global mname_ref
    global nlay, nrow, ncol
    global delr, delc
    global shift_some_x
    global h_start
    global k11
    global chd_spd
    global tops
    global rewet_record
    global wetdry

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
        rewet_record=rewet_record,
        icelltype=1,
        k=hk,
        wetdry=wetdry,
    )

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_ref}.hds",
        budget_filerecord=f"{mname_ref}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_leftmodel(sim):
    global mname_left
    global nlay, nrow, ncol_left
    global delr, delc
    global tops
    global h_start
    global h_left
    global left_chd
    global k11
    global chd_spd_left
    global rewet_record
    global wetdry

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
        save_flows=True,
        rewet_record=rewet_record,
        icelltype=1,
        k=hk,
        wetdry=wetdry,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_left)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_left}.hds",
        budget_filerecord=f"{mname_left}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_rightmodel(sim):
    global mname_right
    global nlay, nrow, ncol_right
    global h_right
    global delr, delc
    global tops
    global h_start
    global right_chd
    global k11
    global shift_x, shift_y
    global chd_spd_right
    global rewet_record
    global wetdry

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
        save_flows=True,
        rewet_record=rewet_record,
        icelltype=1,
        k=hk,
        wetdry=wetdry,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd_right)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname_right}.hds",
        budget_filerecord=f"{mname_right}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_gwfexchange(sim):
    global mname_left, mname_right
    global nrow
    global delc, delr
    global ncol_left

    angldegx = 0.0
    cdist = delr
    gwfgwf_data = [
        [
            (ilay, irow, ncol_left - 1),
            (ilay, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for ilay in range(nlay)
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
        dev_interfacemodel_on=True,
    )


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    print("comparing heads to single model reference...")

    fpth = os.path.join(test.workspace, f"{mname_ref}.hds")
    hds = flopy.utils.HeadFile(fpth)
    fpth = os.path.join(test.workspace, f"{mname_left}.hds")
    hds_l = flopy.utils.HeadFile(fpth)
    fpth = os.path.join(test.workspace, f"{mname_right}.hds")
    hds_r = flopy.utils.HeadFile(fpth)

    times = hds.get_times()
    for iper, t in enumerate(times):
        heads = hds.get_data(totim=t)
        heads_left = hds_l.get_data(totim=t)
        heads_right = hds_r.get_data(totim=t)
        heads_2models = np.append(heads_left, heads_right, axis=2)

        # in this test we want to have the top layer in the left model
        # dry in period 2, but the cells in the right model should remain
        # active. This tests the interface model for dealing with drying
        # and wetting, and handling inactive cells, explicitly
        if iper == 1:
            assert np.all(heads_left[0, 0, :] == -1.0e30), (
                "left model, top layer should be DRY in period 2"
            )
            assert np.all(heads_right[0, 0, :] > -1.0e30), (
                "right model, top layer should be WET in period 2"
            )

        # compare heads
        maxdiff = np.amax(abs(heads - heads_2models))
        assert maxdiff < 10 * hclose_check, (
            f"Max. head diff. {maxdiff} should \
                         be within solver tolerance (x10): {10 * hclose_check}"
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
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

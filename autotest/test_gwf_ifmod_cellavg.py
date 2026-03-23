"""
General test for the interface model approach.
It compares the result of a single reference model
to the equivalent case where the domain is decomposed
and joined by a GWF-GWF exchange.

In this case we test cell averaging options in the exchange,
which is consequently enabled in the interface model such
that the coupled setup leads to identical results.

                  'refmodel'            'leftmodel'       'rightmodel'
               (k = 1.0, k = 0.01)       (k = 1.0)         (k = 0.01)

   layer 1:    1 . . . . . . . 1    =    1 . . . 1    +    1 . . . 1

Hydraulic conductivities are uniform in the sub-models but discontinuous
across the interface in order to test the averaging methods.

We assert (in)equality on the head values. All models are part of the same
solution for convenience. Finally, the budget error is checked.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["avg_hrm", "avg_log", "avg_amt-lmk", "avg_fails_1", "avg_fails_2"]
alt_cellavg_model = ["default", "logarithmic", "amt-lmk", "amt-lmk", "default"]
cellavg_exg = ["harmonic", "logarithmic", "amt-lmk", "harmonic", "logarithmic"]
results_match = [True, True, True, False, False]

# some global convenience...:
# model names
mname_ref = "refmodel"
mname_left = "leftmodel"
mname_right = "rightmodel"

# solver criterion
hclose_check = 1e-9
max_inner_it = 300
nper = 1

# model spatial discretization
nlay = 1
nrow = 1
ncol = 15
ncol_left = 10
ncol_right = 5

delr = 1.0
delc = 1.0
area = delr * delc

# shift (hor. and vert.)
shift_x = ncol_left * delr
shift_y = 0.0

# top/bot of the aquifer
tops = [10.0, 0.0]

# hydraulic conductivity
hk = np.ones((nlay, nrow, ncol), dtype=float)
hk[:, :, ncol_left:] = 0.01
hk_left = hk[:, :, :ncol_left]
hk_right = hk[:, :, ncol_left:]

# boundary stress period data
h_left = 10.0
h_right = 0.0

# initial head
h_start = 0.0

# head boundaries
lchd = [[(0, 0, 0), h_left]]
rchd = [[(0, 0, ncol - 1), h_right]]
rchd_right = [[(0, 0, ncol_right - 1), h_right]]

chd_ref = lchd + rchd
chd_spd = {0: chd_ref}
chd_spd_left = {0: lchd}
chd_spd_right = {0: rchd_right}


def get_model(idx, ws):
    name = cases[idx]
    alt_avg_model = alt_cellavg_model[idx]
    avg_exg = cellavg_exg[idx]

    # parameters and spd
    # tdis
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, max_inner_it
    hclose, rclose, relax = hclose_check, 1e-3, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        inner_dvclose=hclose,
        outer_dvclose=hclose,
    )

    # the full gwf model as a reference
    gwf = add_refmodel(sim)
    if alt_avg_model != "default":
        gwf.npf.alternative_cell_averaging = alt_avg_model

    # now add two coupled models with the interface model enabled,
    # to be stored in the same solution as the reference model
    gwf_left = add_leftmodel(sim)
    if alt_avg_model != "default":
        gwf_left.npf.alternative_cell_averaging = alt_avg_model
    gwf_right = add_rightmodel(sim)
    if alt_avg_model != "default":
        gwf_right.npf.alternative_cell_averaging = alt_avg_model
    gwf_exg = add_gwfexchange(sim)
    gwf_exg.cell_averaging = avg_exg

    return sim


def add_refmodel(sim):
    global mname_ref
    global nlay, nrow, ncol
    global delr, delc
    global h_start
    global chd_spd
    global tops

    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname_ref, save_flows=True)

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

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        k=hk,
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
    global chd_spd_left

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
        k=hk_left,
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
    global shift_x, shift_y
    global chd_spd_right

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
        k=hk_right,
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

    return gwfgwf


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    match_str = "" if results_match[idx] else "not "
    print(
        f"comparing heads to single model reference for "
        f"{alt_cellavg_model[idx]} and {cellavg_exg[idx]}. "
        f"Results should {match_str}match."
    )

    fpth = os.path.join(test.workspace, f"{mname_ref}.hds")
    hds = flopy.utils.HeadFile(fpth)
    fpth = os.path.join(test.workspace, f"{mname_left}.hds")
    hds_l = flopy.utils.HeadFile(fpth)
    fpth = os.path.join(test.workspace, f"{mname_right}.hds")
    hds_r = flopy.utils.HeadFile(fpth)

    t = hds.get_times()[-1]

    heads = hds.get_data(totim=t)
    heads_left = hds_l.get_data(totim=t)
    heads_right = hds_r.get_data(totim=t)
    heads_2models = np.append(heads_left, heads_right, axis=2)

    # compare heads
    maxdiff = np.amax(abs(heads - heads_2models))
    if results_match[idx]:
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
                        f"Cumulative balance error = {cumul_balance_error} "
                        f"for {mname}, should equal 0.0"
                    )
    else:
        assert maxdiff > 10 * hclose_check, (
            f"A significant difference in head should be observed here, "
            f"but it is only {maxdiff}."
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

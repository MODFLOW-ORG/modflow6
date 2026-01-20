"""
General test for the interface model approach.
It compares the skewed decomposition of the domain
to the trivial analytical result (constant gradient).
In this case with the use of idomain to deactivate half
of the sub-models. Note that the cells with idomain==0
overlap with the active cells of the other model.

     'leftmodel'         'rightmodel'

    1 1 1 0 0 0 0       1 1 1 1 1 1 1
    1 1 1 1 0 0 0       0 1 1 1 1 1 1
    1 1 1 1 1 0 0   +   0 0 1 1 1 1 1
    1 1 1 1 1 1 0       0 0 0 1 1 1 1
    1 1 1 1 1 1 1       0 0 0 0 1 1 1

We assert equality on the head values and check budgets.
"""

import os

import flopy
import numpy as np
import pytest
from flopy.mf6.utils import Mf6Splitter
from framework import TestFramework

cases = ["ifmod_skewed"]

# some global convenience...:
# model name
mname = "skewed"

# solver criterion
hclose_check = 1e-9
max_inner_it = 300
nper = 1

# model spatial discretization
nlay = 1
ncol = 10
nrow = 5

# idomain
idomain = np.ones((nlay, nrow, ncol))

delr = 1.0
delc = 1.0
area = delr * delc

# top/bot of the aquifer
tops = [1.0, 0.0]

# hydraulic conductivity
hk = 10.0

# boundary stress period data
h_left = 10.0
h_right = 1.0

# initial head
h_start = 0.0

# head boundaries
lchd = [[(ilay, irow, 0), h_left] for irow in range(nrow) for ilay in range(nlay)]
rchd = [
    [(ilay, irow, ncol - 1), h_right] for irow in range(nrow) for ilay in range(nlay)
]
chd = lchd + rchd

chd_spd = {0: chd}


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

    gwf = flopy.mf6.ModflowGwf(sim, modelname=mname, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        xorigin=0.0,
        yorigin=0.0,
        top=tops[0],
        botm=tops[1:],
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{mname}.hds",
        budget_filerecord=f"{mname}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # split the model
    splitter = Mf6Splitter(sim)
    mask = np.zeros(shape=(nrow, ncol))
    for irow in range(nrow):
        istart = irow + 3
        mask[irow, istart:] = 1
    split_sim = splitter.split_model(mask)
    split_sim.set_sim_path(dir)

    return split_sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    print("comparing heads to single model reference...")

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)

    mname_left = sim.model_names[0]
    mname_right = sim.model_names[1]

    fpth = os.path.join(test.workspace, f"{mname_left}.hds")
    hds_left = flopy.utils.HeadFile(fpth).get_alldata()
    hds_left[hds_left == 1.0e30] = 0.0

    fpth = os.path.join(test.workspace, f"{mname_right}.hds")
    hds_right = flopy.utils.HeadFile(fpth).get_alldata()
    hds_right[hds_right == 1.0e30] = 0.0

    hds = np.zeros((nrow, ncol), dtype=float)
    hds[:, 0:7] = hds[:, 0:7] + hds_left[:, :]
    hds[:, 3:] = hds[:, 3:] + hds_right[:, :]

    cst_gradient = np.linspace(10.0, 1.0, ncol)
    for irow in range(nrow):
        assert hds[irow, :] == pytest.approx(cst_gradient, rel=10 * hclose_check), (
            f"Head values for row {irow} do not match analytical result. "
            f"Expected {cst_gradient}, but got {hds[irow, :]}"
        )

    # check budget error from .lst file
    for mname in [mname_left, mname_right]:
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

"""
Test the specific discharge calculation

  1. for a DIS grid with a chain of cells extending in the y direction
  2. that same grid split in two

i.e.

    [1]         [1]
    ...         ...
    [5]         [5]
           vs    +
    [6]         [1]
    ...         ...
    [10         [5]

Case 2. comes with a non-zero angledegx in the exchange file (=270.0)
and, as a result from rounding errors, a normal vector of the
form (epsilon,1.0) with epsilon tiny (~1.e-16) but not zero (as
it would be for the internal connection of the DIS grid).
The calculation of specific discharge should be able to deal with
these roundoff errors.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

from flopy.mf6.utils import Mf6Splitter

single_case_name = "dis_single"
split_case_name = "dis_split"
cases = [single_case_name, split_case_name]

model_name = "gwf_model"

# solver criterion
hclose = 1e-9
max_inner_it = 50
nper = 1

# model spatial discretization
nlay = 1
nrow = 10
ncol = 1

# idomain
idomain = np.ones((nlay, nrow, ncol))

# cell size
delr = 1.0
delc = 1.0
area = delr * delc

# top/bot of the aquifer
tops = [10.0, 0.0]

# hydraulic conductivity
hk = 10.0

# boundary stress period data
h_north = 5.0
h_south = 0.0

# initial head
h_start = 0.0


# head boundaries
chd_spd = {0: [[(0, 0, 0), h_north], [(0, nrow - 1, 0), h_south]]}

def get_model(idx, dir):
    name = cases[idx]

    # parameters and spd
    # tdis
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))
        
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=50,
        inner_dvclose=0.1*hclose,
        rcloserecord=0.001,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=0.0,
        filename="gwf.ims",
    )

    gwf = flopy.mf6.ModflowGwf(sim, modelname=model_name, save_flows=True)

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
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
        head_filerecord=f"{model_name}.hds",
        budget_filerecord=f"{model_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    if cases[idx] == split_case_name:
        # split
        splitter = Mf6Splitter(sim)
        mask = np.zeros(shape=(nrow, ncol))
        mask[nrow//2:,:] = 1
        split_sim = splitter.split_model(mask)
        split_sim.set_sim_path(dir)

        # gwfgwf = split_sim.get_package("gwfgwf")
        # gwfgwf.dev_interfacemodel_on = True
        sim = split_sim

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)
    return sim, None


def check_output(idx, test):
    print("comparing heads to single model reference...")

    qx = None
    qy = None
    qz = None
    if cases[idx] == single_case_name:
        sim = test.sims[0]
        gwf = sim.get_model(model_name)
        bud = gwf.output.budget()
        spdis = bud.get_data(text="DATA-SPDIS")[0]
        qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    elif cases[idx] == split_case_name:
        sim = test.sims[0]
        gwf0 = sim.get_model(sim.model_names[0])
        gwf1 = sim.get_model(sim.model_names[1])

        bud0 = gwf0.output.budget()
        spdis0 = bud0.get_data(text="DATA-SPDIS")[0]
        qx0, qy0, qz0 = flopy.utils.postprocessing.get_specific_discharge(spdis0, gwf0)

        bud1 = gwf1.output.budget()
        spdis1 = bud1.get_data(text="DATA-SPDIS")[0]
        qx1, qy1, qz1 = flopy.utils.postprocessing.get_specific_discharge(spdis1, gwf1)
        
        qx = np.concatenate((qx0[0,:,0], qx1[0,:,0]))
        qy = np.concatenate((qy0[0,:,0], qy1[0,:,0]))
        qz = np.concatenate((qz0[0,:,0], qz1[0,:,0]))

    # drop first and last nodes (CHD)
    qx = qx[1:-2]
    qy = qy[1:-2]
    qz = qz[1:-2]

    qy_theory = -hk * (h_north - h_south)/ ((nrow - 1) * delc)
    assert np.allclose(qx, 0.0), "spdis cannot have x component in this problem"
    assert np.allclose(qy, qy_theory), "spdis y component should equal theory"
    assert np.allclose(qz, 0.0), "spdis cannot have z component in this problem"

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

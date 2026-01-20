import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "gwf_pertim",
]

# static model data
# temporal discretization
nper = 1
perlen = [0.0]
nstp = [1]
tsmult = [1.0]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# spatial discretization data
nlay, nrow, ncol = 3, 21, 20
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 500.0, 500.0
top = 330.0
botm = [220.0, 200.0, 0.0]
strt = 330.0

# calculate hk
hk = [50.0, 0.01, 200.0]
k33 = [10.0, 0.01, 20.0]

# chd data
canal_spd = [(0, i, 0, 330.0, "canal") for i in range(nrow)]
river_spd = [(0, i, ncol - 1, 320.0, "river") for i in range(nrow)]


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    ims = flopy.mf6.ModflowIms(sim, print_option="ALL", complexity="simple")
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

    # create iterative model solution and register the gwf model with it

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        k=hk,
        k33=k33,
    )

    # chd files
    chd_canal = flopy.mf6.ModflowGwfchd(
        gwf,
        boundnames=True,
        stress_period_data=canal_spd,
        pname="CHD-CANAL",
        filename=f"{name}_canal.chd",
    )
    chd_river = flopy.mf6.ModflowGwfchd(
        gwf,
        boundnames=True,
        stress_period_data=river_spd,
        pname="CHD-RIVER",
        filename=f"{name}_river.chd",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[
            (("BUDGET", "ALL")),
        ],
    )

    return sim, None


def check_output(idx, test):
    fpth = os.path.join(test.workspace, f"{test.name}.lst")
    mflist = flopy.utils.Mf6ListBudget(fpth)
    inc = mflist.get_incremental()

    q_in = 99928.4941
    q_out = 99928.5036
    q_in_sim = inc["CHD_IN"]
    q_out_sim = inc["CHD2_OUT"]

    assert np.allclose([q_in_sim], [q_in]), f"CHD_IN <> {q_in} ({q_in_sim})"
    assert np.allclose([q_out_sim], [q_out]), f"CHD2_OUT <> {q_out} ({q_out_sim})"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

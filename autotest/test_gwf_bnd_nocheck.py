"""
Test to make sure that mf6 is failing with head-dependent boundary conditions with
elevations below the bottom of the layer (only ghb and drain are tested). Also check
that the no check option works with head-dependent boundary conditions with elevations
less than the cell bottom.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

xfail = [True, True, False, False]
boundaries = ["ghb", "drn", "ghb", "drn"]
check = [None, None, True, True]
cases = [f"{bnd}{idx:02d}" for idx, bnd in enumerate(boundaries)]


def build_models(idx, test):
    ws = test.workspace

    bnd = boundaries[idx]

    nlay, nrow, ncol = 2, 1, 1
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0, -1.0]
    strt = 1
    hk = 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        nocheck=check[idx],
    )

    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        print_input=True,
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="simple",
        outer_dvclose=1e-8,
        outer_maximum=200,
        inner_dvclose=1e-9,
        inner_maximum=100,
    )
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="feet",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=hk,
    )

    spd = [((1, 0, 0), strt)]
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(spd),
        stress_period_data=spd,
    )

    spd = [((0, 0, 0), botm[0] - 0.1, 1.0, "bc")]
    if bnd == "ghb":
        p = flopy.mf6.ModflowGwfghb(
            gwf,
            boundnames=True,
            maxbound=len(spd),
            stress_period_data=spd,
        )
    elif bnd == "drn":
        p = flopy.mf6.ModflowGwfdrn(
            gwf,
            boundnames=True,
            maxbound=len(spd),
            stress_period_data=spd,
        )
    obs_data = [
        ("obs1", f"{bnd}", "bc"),
    ]
    obs_name = "bnd.obs"
    p.obs.initialize(
        filename=obs_name,
        continuous={f"{obs_name}.csv": obs_data},
    )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("BUDGET", "ALL"), ("HEAD", "ALL")],
    )

    return sim, None


def check_results(idx, test):
    if not xfail[idx]:
        ws = test.workspace
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
        gwf = sim.get_model()

        bnd = boundaries[idx]
        if bnd == "ghb":
            obs = gwf.ghb.output.obs().get_data()
        elif bnd == "drn":
            obs = gwf.drn.output.obs().get_data()
        v = obs[0][1]
        print(obs)
        answer = -0.55
        assert np.allclose(v, answer), f"boundary flux ({v}) not close to {answer}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_results(idx, t),
        xfail=xfail[idx],
    )
    test.run()

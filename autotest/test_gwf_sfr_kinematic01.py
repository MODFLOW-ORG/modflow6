"""
This test sfr kinematic wave approximation using example
in Ponce (1989) Example 9-3 part 1. In Ponce (1989) the
Courant number was specified, was invariant, and was used
directly in the kinematic wave equation. As a result, the
of the sfr result is approximate since the sfr kinematic
wave approximation is a numerical method where the Courant
number is a function of the solution.

The test has 1 sfr reach that is not connected to the
groundwater model.

The test evaluates EXT-OUTFLOW for reaches 1 against
the known numerical solution for storage weight values
set to 0.5 and 1.0.

Ponce, V. M. (1989). Engineering Hydrology, Principles and Practices.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

paktest = "sfr"
cases = ("sfr-kwe01", "sfr-kwe02")
storage_weights = (0.5, 1.0)


def build_models(idx, test):
    name = cases[idx]

    dt = 60.0 * 60.0
    flows = np.array(
        [
            0.0,
            30.0,
            60.0,
            90.0,
            120.0,
            150.0,
            120.0,
            90.0,
            60.0,
            30.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    times = np.arange(0.0, (flows.shape[0]) * dt, dt)

    # static model data
    # temporal discretization
    nper = times.shape[0] - 1
    nstp = 1
    perioddata = [(dt, nstp, 1.0) for idx in range(nper)]

    # spatial discretization data
    nlay, nrow, ncol = 1, 1, 6
    delr, delc = 100.0, 100.0
    dx = 7200.0
    top = 100.0
    botm = 0.0

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=test.workspace,
        version="mf6",
        exe_name="mf6",
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="seconds", nper=nper, perioddata=perioddata
    )
    ims = flopy.mf6.ModflowIms(sim)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="meters",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=dx,
        delc=dx,
        top=top,
        botm=botm,
    )
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    ic = flopy.mf6.ModflowGwfic(gwf, strt=top)
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-6, sy=0.2, transient={0: True})
    ghb = flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[(0, 0, 0, top, 5.0)])
    oc = flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])

    # sfr file
    nreaches = 1
    slope = 1.0 / dx
    roughness = 0.03574737676661647

    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #        <man> <ncon> <ustrf> <ndv> [<aux(naux)>] [<boundname>]
    pak_data = [
        (ifno, -1, -1, -1, dx, 10.0, slope, top, 1.0, 0.0, roughness, 0, 0.0, 0)
        for ifno in range(nreaches)
    ]
    sfr_spd = {idx: [(0, "inflow", q)] for idx, q in enumerate(flows[1:])}
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_flows=True,
        storage=True,
        dev_storage_weight=storage_weights[idx],
        nreaches=nreaches,
        packagedata=pak_data,
        connectiondata=[(0,)],
        perioddata=sfr_spd,
        pname="sfr-1",
    )
    fname = f"{name}.sfr.obs"
    sfr_obs = {
        f"{fname}.csv": [
            ("inflow", "ext-inflow", (0,)),
            ("outflow", "ext-outflow", (0,)),
        ]
    }
    sfr.obs.initialize(filename=fname, print_input=True, continuous=sfr_obs)

    return sim


def check_output(idx, test):
    print("Checking sfr external outflow")
    name = cases[idx]
    obs_values = flopy.utils.Mf6Obs(test.workspace / f"{name}.sfr.obs.csv").get_data()
    # fmt: off
    test_values = {storage_weights[0]: np.array(
        [ 
         0.        ,  -10.68920803,  -58.58988171,  -92.59647491, -125.23616408, 
         -144.76383592, -117.40352509,  -89.57181956, -64.81693279,  -53.46056541,  
         -11.22910318,   -4.95517783, -2.72890377
         ]
        ),
        storage_weights[1]: np.array(
        [ 
         0.        ,  -15.3918031 ,  -58.58988171,  -92.59647491, -125.23616408, 
         -144.76383592, -117.40352509,  -88.84252125, -62.90365264,  -47.01932528,  
         -19.87588927,   -9.96585666, -5.63155787
         ]
        ),
    }
    assert np.allclose(
            obs_values["OUTFLOW"], test_values[storage_weights[idx]]
        ), f"failed comparison for '{name}' observation"


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

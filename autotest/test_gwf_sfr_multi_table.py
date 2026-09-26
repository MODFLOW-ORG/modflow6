"""Two connected SFR reaches with different static CROSSSECTIONS tables
(reach 0: 2-column, no MANFRACTION; reach 1: 3-column, non-uniform
MANFRACTION), verifying per-reach table matching. Since flow is
conserved between the reaches, expected DEPTH at each is computed
independently via Manning's equation (cross_section_functions.get_depths).
"""

import os

import flopy
import numpy as np
import pytest
from cross_section_functions import get_depths
from framework import TestFramework

cases = ["sfr-multi-table"]

nlay, nrow, ncol = 1, 1, 1
top = 0.0
botm = -10.0
strt = 0.0

nreaches = 2
rlen = 50.0
rwid = 10.0
rough = 0.03
slope = 0.001
rbth = 1.0
rhk = 0.0
ustrf = 1.0
ndv = 0
inflow = 500.0

# reach 0: 3-point V-shape, no MANFRACTION (2-column table)
r0_x = np.array([0.0, 0.5, 1.0]) * rwid
r0_h = np.array([1.0, 0.0, 1.0])
r0_n = np.full(r0_x.shape, rough)

# reach 1: 5-point W-shape, non-uniform MANFRACTION (3-column table)
r1_x = np.array([0.0, 0.2, 0.5, 0.7, 1.0]) * rwid
r1_h = np.array([1.0, 0.0, 0.5, 0.0, 1.0])
r1_manfraction = np.array([1.0, 1.0, 1.5, 1.5, 1.0])
r1_n = rough * r1_manfraction


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 5, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY", complexity="simple")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), strt]])

    packagedata = [
        [0, "none", rlen, rwid, slope, top, rbth, rhk, rough, 1, ustrf, ndv],
        [1, "none", rlen, rwid, slope, top, rbth, rhk, rough, 1, ustrf, ndv],
    ]
    connectiondata = [
        [0, -1],
        [1, 0],
    ]

    tab0_filename = f"{name}.r0.tab"
    tab1_filename = f"{name}.r1.tab"

    budpth = f"{name}.sfr.cbc"
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        nreaches=nreaches,
        packagedata=packagedata,
        crosssections=[[0, tab0_filename], [1, tab1_filename]],
        connectiondata=connectiondata,
        perioddata={0: [[0, "inflow", inflow]]},
        budget_filerecord=budpth,
        pname="SFR-1",
    )
    flopy.mf6.ModflowUtlsfrtab(
        gwf,
        nrow=r0_x.shape[0],
        ncol=2,
        table=[[x, h] for x, h in zip(r0_x / rwid, r0_h)],
        filename=tab0_filename,
        pname="R0TAB",
        parent_file=sfr,
    )
    flopy.mf6.ModflowUtlsfrtab(
        gwf,
        nrow=r1_x.shape[0],
        ncol=3,
        table=[[x, h, m] for x, h, m in zip(r1_x / rwid, r1_h, r1_manfraction)],
        filename=tab1_filename,
        pname="R1TAB",
        parent_file=sfr,
    )

    obs_file = f"{name}.sfr.obs"
    csv_file = obs_file + ".csv"
    obs_dict = {
        csv_file: [
            ("depth0", "depth", (0,)),
            ("depth1", "depth", (1,)),
            ("outflow", "ext-outflow", (1,)),
        ]
    }
    sfr.obs.initialize(
        filename=obs_file, digits=10, print_input=True, continuous=obs_dict
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        saverecord=[("BUDGET", "ALL")],
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    obs_pth = os.path.join(test.workspace, f"{name}.sfr.obs.csv")
    obs = flopy.utils.Mf6Obs(obs_pth).get_data()

    q = np.abs(obs["OUTFLOW"][-1])

    expected_d0 = get_depths(q, r0_x, r0_h, roughness=r0_n, slope=slope)[0]
    expected_d1 = get_depths(q, r1_x, r1_h, roughness=r1_n, slope=slope)[0]
    assert not np.isclose(expected_d0, expected_d1), (
        "test's two cross-sections must produce distinguishably different "
        "depths for this test to be meaningful"
    )

    d0 = obs["DEPTH0"][-1]
    d1 = obs["DEPTH1"][-1]

    assert np.isclose(d0, expected_d0), (
        f"reach 0 depth {d0} != expected {expected_d0} "
        "(wrong table applied to reach 0?)"
    )
    assert np.isclose(d1, expected_d1), (
        f"reach 1 depth {d1} != expected {expected_d1} "
        "(wrong table applied to reach 1, or MANFRACTION dropped?)"
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

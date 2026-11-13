"""
Test advanced transport observations that give a feature number (ID) or a
connection number (ID2). ID must be between 1 and the number of features,
and a feature must have the connection observed. Otherwise the simulation
stops with an error.

Reach 1 is connected to the aquifer, reach 2 is not.
"""

import flopy
import pytest
from framework import TestFramework

cases = {
    "no-conn": (
        [("x", "sft", 2, 1)],
        "specifies iconn = 1, but this is not a valid connection for ID 2.",
    ),
    "id0-sft": ([("x", "sft", 0, 1)], "was assigned ID = 0. ID must be >= 1 and <= 2."),
    "id0-conc": (
        [("x", "concentration", 0)],
        "was assigned ID = 0. ID must be >= 1 and <= 2.",
    ),
    "valid": (
        [("c1", "concentration", 1), ("c2", "concentration", 2), ("q", "sft", 1, 1)],
        None,
    ),
}


def build_models(name, test):
    obs, _ = cases[name]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)

    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    ims_gwf = flopy.mf6.ModflowIms(sim, filename="gwf.ims")
    sim.register_ims_package(ims_gwf, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 2), 8.0]])
    flopy.mf6.ModflowGwfsfr(
        gwf,
        pname="sfr-1",
        save_flows=True,
        nreaches=2,
        packagedata=[
            [0, (0, 0, 0), 10.0, 1.0, 0.001, 9.0, 1.0, 1.0, 0.03, 1, 1.0, 0],
            [1, "none", 10.0, 1.0, 0.001, 9.0, 1.0, 1.0, 0.03, 1, 1.0, 0],
        ],
        connectiondata=[[0, -1], [1, 0]],
        perioddata=[[0, "inflow", 10.0]],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord="gwf.cbc")

    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    ims_gwt = flopy.mf6.ModflowIms(
        sim, filename="gwt.ims", linear_acceleration="bicgstab"
    )
    sim.register_ims_package(ims_gwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt, nlay=1, nrow=1, ncol=3, delr=10.0, delc=10.0, top=10.0, botm=0.0
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtsft(
        gwt,
        flow_package_name="sfr-1",
        packagedata=[[0, 1.0], [1, 1.0]],
        observations={"obs.csv": obs},
    )
    flopy.mf6.ModflowGwtoc(gwt)

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )
    return sim, None


def check_output(name, test):
    _, error = cases[name]
    if error is None:
        assert (test.workspace / "obs.csv").is_file()
    else:
        output = " ".join(" ".join(test.buffs[0]).split())
        assert error in output


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(name, t),
        check=lambda t: check_output(name, t),
        xfail=cases[name][1] is not None,
    )
    test.run()

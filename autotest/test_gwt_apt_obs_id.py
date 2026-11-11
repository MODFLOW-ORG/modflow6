"""
Test advanced transport observations with an invalid feature or connection.

An observation of a feature's connections (LKT, SFT, MWT, UZT and the energy
equivalents) on a feature with no rows in that budget term used to read an
uninitialized index. An ID of 0 also passed the range check. Both are now
errors.

Lake 1 has one connection to the aquifer, lake 2 has none.
"""

import flopy
import pytest


def build(ws, exe, obs):
    sim = flopy.mf6.MFSimulation(sim_ws=str(ws), exe_name=exe)
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
    flopy.mf6.ModflowGwflak(
        gwf,
        pname="lak-1",
        save_flows=True,
        nlakes=2,
        packagedata=[[0, 9.0, 1], [1, 9.0, 0]],
        connectiondata=[
            [0, 0, (0, 0, 0), "horizontal", 0.1, 5.0, 10.0, 5.0, 10.0],
        ],
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
    flopy.mf6.ModflowGwtlkt(
        gwt,
        flow_package_name="lak-1",
        packagedata=[[0, 1.0], [1, 1.0]],
        observations={"obs.csv": obs},
    )
    flopy.mf6.ModflowGwtoc(gwt)

    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )
    sim.write_simulation(silent=True)
    return sim


def run(sim):
    success, buff = sim.run_simulation(silent=True, report=True)
    return success, " ".join(" ".join(buff).split())


def test_obs_feature_without_connections(function_tmpdir, targets):
    sim = build(function_tmpdir, targets["mf6"], [("x", "lkt", 2, 1)])
    success, output = run(sim)
    assert not success
    assert "specifies iconn = 1, but this is not a valid connection for ID 2." in (
        output
    )


@pytest.mark.parametrize("obs", [("x", "lkt", 0, 1), ("x", "concentration", 0)])
def test_obs_id_zero(function_tmpdir, targets, obs):
    sim = build(function_tmpdir, targets["mf6"], [obs])
    success, output = run(sim)
    assert not success
    assert "was assigned ID = 0. ID must be >= 1 and <= 2." in output


def test_obs_valid(function_tmpdir, targets):
    obs = [("c1", "concentration", 1), ("c2", "concentration", 2), ("q", "lkt", 1, 1)]
    sim = build(function_tmpdir, targets["mf6"], obs)
    success, output = run(sim)
    assert success, output

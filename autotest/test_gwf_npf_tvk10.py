"""
Test that a time-varying hydraulic conductivity (TVK) value given in a
PERIOD block at an inactive cell (IDOMAIN 0) is an input error.
"""

import flopy

name = "tvk10"
inactive = (0, 0, 2)


def test_inactive_cell(function_tmpdir, targets):
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=function_tmpdir, exe_name=targets["mf6"]
    )
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=4, top=10.0, botm=0.0, idomain=[1, 1, 0, 1]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowUtltvk(
        npf, perioddata={1: [((0, 0, 1), "K", 2.0)]}, filename=f"{name}.tvk"
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 10.0], [(0, 0, 3), 5.0]]
    )
    sim.write_simulation(silent=True)

    # flopy rejects an inactive cell, so write the PERIOD block directly
    lay, row, col = (i + 1 for i in inactive)
    (function_tmpdir / f"{name}.tvk").write_text(
        f"BEGIN period 2\n  {lay} {row} {col} K 2.0\nEND period 2\n"
    )

    success, buff = sim.run_simulation(silent=True, report=True)
    assert not success, "simulation with a TVK value at an inactive cell ran"
    text = "\n".join(buff)
    assert "is not in the active model domain" in text, text

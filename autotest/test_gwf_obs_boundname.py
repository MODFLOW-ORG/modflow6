"""
Test that an observation can't name a boundary in a package without BOUNDNAMES.

The name array has no elements then, so the search read out of bounds, and a list
package wrote the no-data value with no error.
"""

import flopy
import pytest

errmsg = "but BOUNDNAMES is not specified"


@pytest.mark.parametrize("boundnames", [True, False])
def test_obs_boundname(function_tmpdir, targets, boundnames):
    sim = flopy.mf6.MFSimulation(sim_ws=str(function_tmpdir), exe_name=targets["mf6"])
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, ncol=3)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 2), 1.0]])
    row = [(0, 0, 0), -1.0] + (["mywell"] if boundnames else [])
    flopy.mf6.ModflowGwfwel(
        gwf,
        boundnames=boundnames,
        stress_period_data=[row],
        observations={"wel.csv": [("x", "wel", "mywell")]},
    )
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    if boundnames:
        assert success, "\n".join(buff)
    else:
        assert not success
        output = " ".join(" ".join(buff).split())
        assert errmsg in output
        assert "gwf.wel" in output

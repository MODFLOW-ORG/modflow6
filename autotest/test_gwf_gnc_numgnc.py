"""
Test that a ghost-node correction's GNCDATA block has exactly NUMGNC rows.

With too few, the rest of GNC's arrays stay unset, and the node numbers of 0 read
out of bounds. With too many, the extra rows are ignored, which changes the heads.
"""

import flopy
import pytest

# cell 1 next to cell 2, with cell 3 contributing
ROW = [(0, 0, 0), (0, 0, 1), (0, 0, 2), 0.25]
cases = {
    "match": (1, [ROW], None),
    "too-high": (2, [ROW], "rows, but NUMGNC is 2"),
    "too-low": (
        1,
        [ROW, [(0, 0, 1), (0, 0, 2), (0, 0, 0), 0.25]],
        "more rows than NUMGNC",
    ),
}


@pytest.mark.parametrize("case", cases)
def test_gnc_numgnc(function_tmpdir, targets, case):
    numgnc, gncdata, errmsg = cases[case]
    sim = flopy.mf6.MFSimulation(sim_ws=str(function_tmpdir), exe_name=targets["mf6"])
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim, linear_acceleration="BICGSTAB", inner_dvclose=1e-6, outer_dvclose=1e-6
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf")
    flopy.mf6.ModflowGwfdis(gwf, ncol=3)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 1.0], [(0, 0, 2), 0.0]]
    )
    flopy.mf6.ModflowGwfgnc(gwf, numgnc=numgnc, numalphaj=1, gncdata=gncdata)
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    if errmsg is None:
        assert success, "\n".join(buff)
    else:
        assert not success
        assert errmsg in " ".join(" ".join(buff).split())

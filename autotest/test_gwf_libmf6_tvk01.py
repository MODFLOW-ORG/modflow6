"""
Test when time-varying hydraulic conductivity (TVK) values are applied,
through the API. In period 2, TVK gives K at cell A as a value and at cell
B by a time series. As in 6.7.0, a value given in a PERIOD block is applied
in prepare_time_step and a time series value in prepare_solve. A K set
through the API at cell A after prepare_time_step is kept.
"""

import flopy
import numpy as np
from modflow_devtools.markers import requires_pkg

name = "libtvk01"
k0, k1, kts, kapi = 1.0, 3.0, 5.0, 7.0
cell_a, node_a = (0, 0, 1), 1  # zero-based node numbers
cell_b, node_b = (0, 0, 2), 2


def build(ws):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=4, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0)
    tvk = flopy.mf6.ModflowUtltvk(
        npf, perioddata={1: [(cell_a, "K", k1), (cell_b, "K", "kts")]}
    )
    tvk.ts.initialize(
        filename=f"{name}.tvk.ts",
        timeseries=[(0.0, k0), (1.0, kts), (2.0, kts)],
        time_series_namerecord="kts",
        interpolation_methodrecord="stepwise",
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), 10.0], [(0, 0, 3), 5.0]]}
    )
    sim.write_simulation(silent=True)


def solve(mf6):
    """Solve the time step after prepare_solve and finalize it."""
    mxiter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))[0]
    for _ in range(mxiter):
        if mf6.solve(1):
            break
    mf6.finalize_solve(1)
    mf6.finalize_time_step()


@requires_pkg("modflowapi")
def test_tvk_api_timing(function_tmpdir, targets):
    from modflowapi import ModflowApi

    build(function_tmpdir)
    mf6 = ModflowApi(targets["libmf6"], working_directory=function_tmpdir)
    mf6.initialize()
    try:
        k11_addr = mf6.get_var_address("K11", name.upper(), "NPF")

        # period 1
        mf6.prepare_time_step(mf6.get_time_step())
        mf6.prepare_solve(1)
        solve(mf6)

        # period 2: the value at cell A is applied, the series at B is not
        mf6.prepare_time_step(mf6.get_time_step())
        k11 = mf6.get_value(k11_addr)
        assert np.allclose(k11[[node_a, node_b]], [k1, k0]), (
            f"K {k11[[node_a, node_b]]} after prepare_time_step, expected {[k1, k0]}"
        )

        # an API value at cell A is kept; the series value is applied at B
        k11[node_a] = kapi
        mf6.set_value(k11_addr, k11)
        mf6.prepare_solve(1)
        k11 = mf6.get_value(k11_addr)
        assert np.allclose(k11[[node_a, node_b]], [kapi, kts]), (
            f"K {k11[[node_a, node_b]]} after prepare_solve, expected {[kapi, kts]}"
        )
        solve(mf6)
    finally:
        mf6.finalize()

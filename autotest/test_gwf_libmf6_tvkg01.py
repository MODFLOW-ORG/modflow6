"""
Test when time-varying hydraulic conductivity grid array (TVKG) values are
applied, through the API. In period 2, TVKG gives K at cell A as an array
and K33 at cell B by a time-array series. As with TVK, an array given in a
PERIOD block is applied in prepare_time_step and a time-array series value
in prepare_solve. A K set through the API at cell A after
prepare_time_step is kept.
"""

import flopy
import numpy as np
from modflow_devtools.markers import requires_pkg

name = "libtvkg01"
DNODATA = 3.0e30
ncol = 4
k0, k1, kapi = 1.0, 3.0, 7.0
k330, k33ts = 1.0, 0.2
node_a, node_b = 1, 2  # zero-based node numbers


def grid(values):
    """Full grid array, DNODATA except at the given nodes."""
    a = np.full((1, 1, ncol), DNODATA)
    for n, v in values.items():
        a[0, 0, n] = v
    return a


def build(ws):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=ncol, top=10.0, botm=0.0)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=k330)
    npf.tvkg.initialize(
        filename=f"{name}.tvkg",
        k={1: grid({node_a: k1})},
        k33={1: "TIMEARRAYSERIES k33tas"},
    )
    npf.tvkg.tas.initialize(
        filename=f"{name}.tvkg.tas",
        tas_array={t: grid({node_b: k33ts}) for t in (0.0, 1.0, 2.0)},
        time_series_namerecord="k33tas",
        interpolation_methodrecord="stepwise",
    )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), 10.0], [(0, 0, ncol - 1), 5.0]]}
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
def test_tvkg_api_timing(function_tmpdir, targets):
    from modflowapi import ModflowApi

    build(function_tmpdir)
    mf6 = ModflowApi(targets["libmf6"], working_directory=function_tmpdir)
    mf6.initialize()
    try:
        k11_addr = mf6.get_var_address("K11", name.upper(), "NPF")
        k33_addr = mf6.get_var_address("K33", name.upper(), "NPF")

        # period 1
        mf6.prepare_time_step(mf6.get_time_step())
        mf6.prepare_solve(1)
        solve(mf6)

        # period 2: the array at cell A is applied, the series at B is not
        mf6.prepare_time_step(mf6.get_time_step())
        k11, k33 = mf6.get_value(k11_addr), mf6.get_value(k33_addr)
        assert np.isclose(k11[node_a], k1), (
            f"K {k11[node_a]} after prepare_time_step, expected {k1}"
        )
        assert np.isclose(k33[node_b], k330), (
            f"K33 {k33[node_b]} after prepare_time_step, expected {k330}"
        )

        # an API value at cell A is kept; the series value is applied at B
        k11[node_a] = kapi
        mf6.set_value(k11_addr, k11)
        mf6.prepare_solve(1)
        k11, k33 = mf6.get_value(k11_addr), mf6.get_value(k33_addr)
        assert np.isclose(k11[node_a], kapi), (
            f"K {k11[node_a]} after prepare_solve, expected {kapi}"
        )
        assert np.isclose(k33[node_b], k33ts), (
            f"K33 {k33[node_b]} after prepare_solve, expected {k33ts}"
        )
        solve(mf6)
    finally:
        mf6.finalize()

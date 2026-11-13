"""
Test that FACTOR does not scale DNODATA (3.0E+30) in time-array series
arrays. RCHA (READASARRAYS) RECHARGE is linked to a stepwise series whose
arrays are read with FACTOR 2.0 and hold DNODATA at an inactive cell. After
the first time step is prepared, RECHARGE holds the scaled values and the
unscaled DNODATA value.
"""

import flopy
import numpy as np
from modflow_devtools.markers import requires_pkg

name = "tasfac"
factor = 2.0
values = [0.001, 0.002, 3.0e30]


def build(ws):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowIms(sim)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, top=10.0, botm=0.0, idomain=[1, 1, 0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 0), 10.0]]})
    flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0, pname="RCH-1", filename=f"{name}.rcha")
    sim.write_simulation(silent=True)

    # package and series files with a series array FACTOR
    (ws / f"{name}.rcha").write_text(
        "BEGIN options\n  READASARRAYS\n"
        f"  TAS6 FILEIN {name}.rcha.tas\nEND options\n\n"
        "BEGIN period 1\n  RECHARGE TIMEARRAYSERIES rts\nEND period 1\n"
    )
    row = " ".join(f"{v:.8e}" for v in values)
    records = "".join(
        f"BEGIN time {t}\n  INTERNAL FACTOR {factor}\n    {row}\nEND time\n\n"
        for t in (0.0, 1.0)
    )
    (ws / f"{name}.rcha.tas").write_text(
        f"BEGIN attributes\n  NAME rts\n  METHOD STEPWISE\nEND attributes\n\n{records}"
    )


@requires_pkg("modflowapi")
def test_tas_factor_dnodata(function_tmpdir, targets):
    from modflowapi import ModflowApi

    build(function_tmpdir)
    mf6 = ModflowApi(targets["libmf6"], working_directory=function_tmpdir)
    mf6.initialize()
    try:
        mf6.prepare_time_step(mf6.get_time_step())
        recharge = mf6.get_value(mf6.get_var_address("RECHARGE", name.upper(), "RCH-1"))
    finally:
        mf6.finalize()
    expected = [v if v == 3.0e30 else v * factor for v in values]
    assert np.allclose(recharge, expected), f"RECHARGE {recharge}, expected {expected}"

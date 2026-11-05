"""
Test of the saved interface elevation (ZETA FILEOUT) of the SWI Package.

The solve uses the interface elevation computed from the heads without
limits. The saved zeta is limited to the saturated part of the cell: the cell
bottom, and the top of the saturated column. That top is the cell top in a
confined cell, and the water table (the freshwater head) in a convertible cell
where the head is below the cell top. Where zeta is above the top of the
saturated column the cell holds no freshwater and saltwater fills the column,
so the interface is reported at the water table, not at the cell top.

A one-row, single-fluid model with three constant-head cells covers the three
outcomes (alphaf = 40, saltwater head 0, top 10 m, bottom -40 m):

  head -0.5 m: zeta = 20 m, above the cell top
  head  0.2 m: zeta = -8 m, inside the cell
  head  2.0 m: zeta = -80 m, below the cell bottom

Cases:
  zeta-convertible : ICELLTYPE 1, first cell saved at the water table (-0.5 m)
  zeta-confined    : ICELLTYPE 0, first cell saved at the cell top (10 m)
  zeta-negative    : ICELLTYPE -1 without THICKSTRT (convertible)
  zeta-thickstrt   : ICELLTYPE -1 with THICKSTRT (confined, starting
                     thickness), first cell saved at the cell top
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["zeta-convertible", "zeta-confined", "zeta-negative", "zeta-thickstrt"]
icelltype = [1, 0, -1, -1]
thickstrt = [False, False, False, True]

top, bot = 10.0, -40.0
heads = [-0.5, 0.2, 2.0]
expected = {
    "zeta-convertible": [-0.5, -8.0, bot],
    "zeta-confined": [top, -8.0, bot],
    "zeta-negative": [-0.5, -8.0, bot],
    "zeta-thickstrt": [top, -8.0, bot],
}


def build_models(idx, test):
    ws = test.workspace
    name = "zeta"
    ncol = len(heads)
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=100.0, delc=1.0, top=top, botm=bot
    )
    flopy.mf6.ModflowGwfic(gwf, strt=heads)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=icelltype[idx], thickstrt=thickstrt[idx], k=10.0
    )
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=0.2, steady_state={0: True})
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=f"{name}.zta")
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, j), h] for j, h in enumerate(heads)]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    return sim, None


def check_output(idx, test):
    ws = test.workspace
    head = flopy.utils.HeadFile(ws / "zeta.hds").get_data().ravel()
    assert np.allclose(head, heads), f"heads {head} differ from the CHD heads"
    zeta = flopy.utils.HeadFile(ws / "zeta.zta", text="zeta").get_data().ravel()
    zexp = expected[cases[idx]]
    assert np.allclose(zeta, zexp), f"{cases[idx]}: saved zeta {zeta}, expected {zexp}"


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

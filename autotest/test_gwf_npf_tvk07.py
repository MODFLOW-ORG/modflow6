"""
Test that a time-varying hydraulic conductivity (TVK) value keeps the value
it was given until TVK changes it again, when K33 is given as a ratio of K
(NPF K33OVERK). A later change of K alone must not change K33.

Steady flow in two layers of a row of cells between constant heads at the
two ends of layer 1. Each case is compared with the same TVK changes and
K33 given as its value (the ratio times K when the ratio is given), without
K33OVERK. Period 3 has two time steps.

Case 1 (tvk07a): K33 ratio at a cell in period 2, K at the same cell in
    period 3.
Case 2 (tvk07b): K33 ratio at a cell in period 2, K at the same cell from a
    time series in period 3 that changes at the second time step.
Case 3 (tvk07c): K from a time series and a K33 ratio at a cell in the same
    period 3 block. The ratio uses K at the first time step of the period
    and is not applied again when the series changes K at the second.
Case 4 (tvk07d): K from a time series at a cell in period 2, replaced by a
    value in period 3. The series changes at the second time step of
    period 3, which must have no effect. This case guards link removal;
    earlier versions also pass it.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk07a", "tvk07b", "tvk07c", "tvk07d"]

nlay, nrow, ncol = 2, 1, 5
k0, k33ratio0 = 1.0, 0.5
cell = (0, 0, 1)
k33ratio1 = 0.2
k1 = 3.0
# time series (time, K); period 3 runs from 2.0 to 3.0 in two steps
kts = {
    "tvk07b": [(0.0, k0), (2.0, k1), (2.5, 5.0), (3.0, 5.0)],
    "tvk07c": [(0.0, k0), (2.0, 2.0), (2.5, 4.0), (3.0, 4.0)],
    "tvk07d": [(0.0, 2.0), (2.5, 6.0), (3.0, 6.0)],
}


def perioddata(name, ratio):
    """TVK period data; the reference gives K33 as its value."""
    if name == "tvk07a":
        k33 = k33ratio1 if ratio else k33ratio1 * k0
        return {1: [(cell, "K33", k33)], 2: [(cell, "K", k1)]}
    if name == "tvk07b":
        k33 = k33ratio1 if ratio else k33ratio1 * k0
        return {1: [(cell, "K33", k33)], 2: [(cell, "K", "kts")]}
    if name == "tvk07c":
        # K at the first time step of period 3 is the series value at 2.0
        k33 = k33ratio1 if ratio else k33ratio1 * 2.0
        return {2: [(cell, "K", "kts"), (cell, "K33", k33)]}
    # tvk07d: the reference gives the series value of period 2 as a value
    if ratio:
        return {1: [(cell, "K", "kts")], 2: [(cell, "K", k1)]}
    return {1: [(cell, "K", 2.0)], 2: [(cell, "K", k1)]}


def build_sim(ws, name, exe, ratio):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim, nper=3, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 2, 1.0)]
    )
    flopy.mf6.ModflowIms(sim, inner_dvclose=1e-10, outer_dvclose=1e-10)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=10.0, botm=[5.0, 0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    if ratio:
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=k33ratio0, k33overk=True)
    else:
        # reference: K33 given as its value
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=k33ratio0 * k0)
    tvk = flopy.mf6.ModflowUtltvk(
        npf, perioddata=perioddata(name, ratio), filename=f"{name}.tvk"
    )
    if name in kts and (ratio or name != "tvk07d"):
        tvk.ts.initialize(
            filename=f"{name}.tvk.ts",
            timeseries=kts[name],
            time_series_namerecord="kts",
            interpolation_methodrecord="stepwise",
        )
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 10.0], [(0, 0, ncol - 1), 5.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    sim.write_simulation(silent=True)
    return sim


def build_models(idx, test):
    return build_sim(test.workspace, cases[idx], "mf6", ratio=True), None


def check_output(idx, test):
    name = cases[idx]

    # reference: K33 given as its value
    ws = test.workspace / "ref"
    sim = build_sim(ws, name, test.targets["mf6"], ratio=False)
    success, _ = sim.run_simulation(silent=True)
    assert success, "reference simulation failed"

    h = flopy.utils.HeadFile(test.workspace / f"{name}.hds").get_alldata()
    href = flopy.utils.HeadFile(ws / f"{name}.hds").get_alldata()
    assert np.abs(href[-1] - href[0]).max() > 1e-3, "TVK changes have no effect"
    for i in range(len(h)):
        assert np.allclose(h[i], href[i], atol=1e-10), (
            f"time step {i + 1} heads differ from the reference: "
            f"max diff {np.abs(h[i] - href[i]).max()}"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

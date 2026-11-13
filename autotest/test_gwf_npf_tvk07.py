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
    period 3 block. A value given in a PERIOD block is applied before a
    time series value, so the ratio uses the K in effect before period 3
    (6.7.0 used the series value of the previous time step, which differs
    here), and is not applied again when the series changes K.
Case 4 (tvk07d): K from a time series at a cell in period 2, replaced by a
    value in period 3. The series changes at the second time step of
    period 3, which must have no effect. This case guards link removal;
    earlier versions also pass it.
Case 5 (tvk07e): a K33 ratio listed before K at the same cell in the
    period 2 block. K is applied first, so the ratio uses the new K
    whatever the order of the lines (6.7.0 used the K from before).
Case 6 (tvk07f): as case 1 with a K22 ratio (NPF K22OVERK). K22 does not
    affect the heads of this one-row model, so the K22 values listed with
    PRINT_INPUT are checked: K22 is set once, from the ratio and K0.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk07a", "tvk07b", "tvk07c", "tvk07d", "tvk07e", "tvk07f"]

nlay, nrow, ncol = 2, 1, 5
k0, k33ratio0 = 2.0, 0.5
cell = (0, 0, 1)
k33ratio1 = 0.2
k22ratio0, k22ratio1 = 0.8, 0.4
k1 = 3.0
# time series (time, K); period 3 runs from 2.0 to 3.0 in two steps
kts = {
    "tvk07b": [(0.0, k0), (2.0, k1), (2.5, 5.0), (3.0, 5.0)],
    "tvk07c": [(0.0, k0), (1.0, 7.0), (2.0, 2.0), (2.5, 4.0), (3.0, 4.0)],
    "tvk07d": [(0.0, 1.5), (2.5, 6.0), (3.0, 6.0)],
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
        # the ratio is applied before the series sets K in period 3
        k33 = k33ratio1 if ratio else k33ratio1 * k0
        return {2: [(cell, "K", "kts"), (cell, "K33", k33)]}
    if name == "tvk07e":
        k33 = k33ratio1 if ratio else k33ratio1 * k1
        return {1: [(cell, "K33", k33), (cell, "K", k1)]}
    if name == "tvk07f":
        k22 = k22ratio1 if ratio else k22ratio1 * k0
        return {1: [(cell, "K22", k22)], 2: [(cell, "K", k1)]}
    # tvk07d: the reference gives the series value of period 2 as a value
    if ratio:
        return {1: [(cell, "K", "kts")], 2: [(cell, "K", k1)]}
    return {1: [(cell, "K", 1.5)], 2: [(cell, "K", k1)]}


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
    if ratio and name == "tvk07f":
        npf = flopy.mf6.ModflowGwfnpf(
            gwf, k=k0, k22=k22ratio0, k22overk=True, k33=k33ratio0, k33overk=True
        )
    elif name == "tvk07f":
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k22=k22ratio0 * k0, k33=k33ratio0 * k0)
    elif ratio:
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=k33ratio0, k33overk=True)
    else:
        # reference: K33 given as its value
        npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=k33ratio0 * k0)
    tvk = flopy.mf6.ModflowUtltvk(
        npf,
        print_input=True,
        perioddata=perioddata(name, ratio),
        filename=f"{name}.tvk",
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

    if name == "tvk07f":
        text = (test.workspace / f"{name}.lst").read_text()
        pattern = r"Setting K22 value for cell .*= +(\S+)"
        k22 = [float(v) for v in re.findall(pattern, text)]
        assert np.allclose(k22, [k22ratio1 * k0]), (
            f"K22 values set {k22}, expected once: {k22ratio1 * k0}"
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

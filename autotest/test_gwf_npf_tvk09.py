"""
Test time-varying hydraulic conductivity (TVK) in two models of one
simulation. The models share the simulation's input step, and each TVK
package applies only its own values. Model A gives a K33 ratio (NPF
K33OVERK) in period 2 and K from a time series in period 3; model B gives
K from a time series in period 2 and K as a value in period 3. The heads of
each model run together must equal the heads of that model run alone.
"""

import flopy
import numpy as np

names = ["tvk09a", "tvk09b"]
nlay, nrow, ncol = 2, 1, 5
k0 = 1.0
cell = (0, 0, 1)
perioddata = {
    "tvk09a": {1: [(cell, "K33", 0.2)], 2: [(cell, "K", "kts")]},
    "tvk09b": {1: [(cell, "K", "kts")], 2: [(cell, "K", 3.0)]},
}
kts = {
    "tvk09a": [(0.0, k0), (2.0, 3.0), (2.5, 5.0), (3.0, 5.0)],
    "tvk09b": [(0.0, k0), (1.0, 2.0), (2.5, 6.0), (3.0, 6.0)],
}


def add_model(sim, name):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, top=10.0, botm=[5.0, 0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=k0, k33=0.5, k33overk=True)
    tvk = flopy.mf6.ModflowUtltvk(
        npf, perioddata=perioddata[name], filename=f"{name}.tvk"
    )
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
    return gwf


def run_sim(ws, models, exe):
    sim = flopy.mf6.MFSimulation(sim_name="tvk09", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim, nper=3, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0), (1.0, 2, 1.0)]
    )
    ims = flopy.mf6.ModflowIms(sim, inner_dvclose=1e-10, outer_dvclose=1e-10)
    gwfs = [add_model(sim, name) for name in models]
    sim.register_ims_package(ims, [gwf.name for gwf in gwfs])
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success, f"simulation of {models} failed"


def test_two_models(function_tmpdir, targets):
    exe = targets["mf6"]
    run_sim(function_tmpdir / "both", names, exe)
    for name in names:
        run_sim(function_tmpdir / name, [name], exe)
        h = flopy.utils.HeadFile(function_tmpdir / "both" / f"{name}.hds")
        href = flopy.utils.HeadFile(function_tmpdir / name / f"{name}.hds")
        h, href = h.get_alldata(), href.get_alldata()
        assert np.abs(href[-1] - href[0]).max() > 1e-3, (
            f"{name}: TVK changes have no effect"
        )
        assert np.allclose(h, href, atol=1e-10), (
            f"{name}: heads with two models differ from one model, "
            f"max diff {np.abs(h - href).max()}"
        )

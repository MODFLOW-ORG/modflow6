"""
Test time-varying hydraulic conductivity (TVK) values when adaptive time
stepping (ATS) retries the first time step. Each retry backs up in time,
so the values written in that step must be applied again. The model is the
one used in test_gwf_ats03, with K at one cell given by a linear time
series. The K value applied for each attempt must be the series value at
the end of that attempt.

Case 1 (tvk08a): K from the series.
Case 2 (tvk08b): K from the series and a K33 ratio (NPF K33OVERK) in the
    same block; each attempt applies the ratio to that attempt's K.

These cases guard the time step retry path; earlier versions that applied
all values at every time step also pass them.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk08a", "tvk08b"]
nlay, nrow, ncol = 1, 1, 10
k0, k1, tend = 1.0, 5.0, 100.0
k33ratio = 0.5


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(100.0, 1, 1.0)]
    )
    # dt0, dtmin, dtmax, dtadj, dtfailadj
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, 100.0, 1.0e-5, 100.0, 2.0, 5.0)],
        filename=f"{name}.ats",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1.5e-6,
        outer_maximum=15,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
        inner_maximum=5,
        inner_dvclose=1.5e-6,
        rcloserecord=1e-6,
        linear_acceleration="CG",
        relaxation_factor=0.97,
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100.0, delc=1.0, top=100.0, botm=0.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=100.0)
    rows = [((0, 0, 4), "K", "kts")]
    if name == "tvk08b":
        npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k0, k33=1.0, k33overk=True)
        rows.append(((0, 0, 4), "K33", k33ratio))
    else:
        npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k0)
    tvk = flopy.mf6.ModflowUtltvk(
        npf, print_input=True, perioddata={0: rows}, filename=f"{name}.tvk"
    )
    tvk.ts.initialize(
        filename=f"{name}.tvk.ts",
        timeseries=[(0.0, k0), (tend, k1)],
        time_series_namerecord="kts",
        interpolation_methodrecord="linearend",
    )
    flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, ss=0.0, sy=0.1, steady_state={0: False}, transient={0: True}
    )
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [[(0, 0, 0), -1.0]]})
    chd = flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, ncol - 1), "tshead"]]}
    )
    chd.ts.initialize(
        filename=f"{name}.chd.ts",
        timeseries=[(0.0, 100.0), (100.0, 50.0)],
        time_series_namerecord="tshead",
        interpolation_methodrecord="linearend",
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    return sim, None


def applied(text, field):
    """Values logged for a field, in order."""
    pattern = rf"Setting {field} value for cell .*= +(\S+)"
    return [float(v) for v in re.findall(pattern, text)]


def check_output(idx, test):
    name = cases[idx]

    # length of each attempt of the first time step
    text = (test.workspace / "mfsim.lst").read_text()
    first = re.findall(r"time step set to\s+(\S+)\s+for step 1 and period 1", text)
    retries = re.findall(
        r"step 1 and period 1 will be retried using time step of\s+(\S+)", text
    )
    assert len(retries) > 0, "no ATS time step was retried"
    dt = [float(v) for v in first[:1] + retries]
    n = len(dt)

    # values applied for each attempt, before the second time step
    text = (test.workspace / f"{name}.lst").read_text()
    kvals = applied(text, "K")
    expected = [k0 + (k1 - k0) * t / tend for t in dt]
    assert np.allclose(kvals[:n], expected), (
        f"K applied for each attempt {kvals[:n]}, expected {expected}"
    )
    if name == "tvk08b":
        k33vals = applied(text, "K33")
        assert np.allclose(k33vals[:n], [k33ratio * k for k in expected]), (
            f"K33 applied for each attempt {k33vals[:n]}, expected "
            f"{[k33ratio * k for k in expected]}"
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

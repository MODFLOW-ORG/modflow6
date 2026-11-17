"""
Test time-varying hydraulic conductivity grid array (TVKG) values with
adaptive time stepping (ATS). The model is the one used in
test_gwf_npf_tvk08, with K at one cell given by a linear time-array series
(DNODATA at the other cells). The first time step is retried; each retry
backs up in time, so the series value is applied again. The K value
applied for each attempt of every time step must be the series average
over that attempt. A second period without ATS follows; its one time step
applies the series value once.

Case 1 (tvkg_atsa): K from the series.
Case 2 (tvkg_atsb): K from the series and a K33 ratio array (NPF
    K33OVERK) in the same block. The ratio array is applied once, before
    the series sets K, so it uses K from before the period, and a retry
    does not apply it again.
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvkg_atsa", "tvkg_atsb"]
nlay, nrow, ncol = 1, 1, 10
perlen = [100.0, 1.0]
k0, k1, tend = 2.0, 5.0, sum(perlen)
k33ratio = 0.5
cell = (0, 0, 4)
DNODATA = 3.0e30


def grid(value):
    """Full grid array, DNODATA except at the cell."""
    a = np.full((nlay, nrow, ncol), DNODATA)
    a[cell] = value
    return a


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(perlen),
        perioddata=[(p, 1, 1.0) for p in perlen],
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
    if name == "tvkg_atsb":
        npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k0, k33=1.0, k33overk=True)
    else:
        npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k0)
    kwargs = {"k": {0: "TIMEARRAYSERIES ktas"}}
    if name == "tvkg_atsb":
        kwargs["k33"] = {0: grid(k33ratio)}
    npf.tvkg.initialize(filename=f"{name}.tvkg", print_input=True, **kwargs)
    npf.tvkg.tas.initialize(
        filename=f"{name}.tvkg.tas",
        tas_array={0.0: grid(k0), tend: grid(k1)},
        time_series_namerecord="ktas",
        interpolation_methodrecord="linear",
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
        timeseries=[(0.0, 100.0), (tend, 100.0 - 0.5 * tend)],
        time_series_namerecord="tshead",
        interpolation_methodrecord="linearend",
    )
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    return sim, None


def applied(text, field):
    """Values TVKG logged for a field, in order."""
    pattern = rf"TVKG: setting {field} for cell .*= +(\S+)"
    return [float(v) for v in re.findall(pattern, text)]


def attempts(text):
    """(start time, length) of each ATS attempt, in order."""
    pattern = r"(time step set to|will be retried using time step of)\s+(\S+)"
    events = [(kind, float(v)) for kind, v in re.findall(pattern, text)]
    result, start = [], 0.0
    for i, (kind, dt) in enumerate(events):
        result.append((start, dt))
        retried = i + 1 < len(events) and events[i + 1][0].startswith("will")
        if not retried:
            start += dt
    return result


def check_output(idx, test):
    name = cases[idx]

    # each attempt of the ATS period, including the retries of step 1
    text = (test.workspace / "mfsim.lst").read_text()
    assert "will be retried" in text, "no ATS time step was retried"
    att = attempts(text)
    t, dt = att[-1]
    assert np.isclose(t + dt, perlen[0]), f"ATS period ended at {t + dt}"

    # value applied for each attempt, then for the step of period 2: the
    # linear series average over it
    att.append((perlen[0], perlen[1]))
    text = (test.workspace / f"{name}.lst").read_text()
    kvals = applied(text, "K")
    assert len(kvals) == len(att), f"K applied {len(kvals)} times, expected {len(att)}"
    expected = [k0 + (k1 - k0) * (t + 0.5 * dt) / tend for t, dt in att]
    # the listing gives five significant digits
    assert np.allclose(kvals, expected, rtol=1.0e-4), (
        f"K applied {kvals}, expected {expected}"
    )
    if name == "tvkg_atsb":
        k33vals = applied(text, "K33")
        assert np.allclose(k33vals, [k33ratio * k0]), (
            f"K33 applied {k33vals}, expected once: {k33ratio * k0}"
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

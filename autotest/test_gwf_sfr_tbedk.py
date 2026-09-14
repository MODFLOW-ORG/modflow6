"""Tests that a PERIOD BEDK/MANNING/UPSTREAM_FRAC override takes
precedence over a reach's own PACKAGEDATA time series, and that the time
series keeps driving the value every timestep until an override is given.

  sfr-bedk           literal BEDK override, reissued every period
  sfr-bedkts         BEDK override given once, persists across later
                      periods without being reissued
  sfr-rhkts          no override; the PACKAGEDATA time series alone
                      drives the value
  sfr-rhkts-nstp     same as sfr-rhkts, with two timesteps per period, to
                      confirm the value updates every timestep and not
                      just once per period
  sfr-manning-ts     same check for MANNING
  sfr-ustrf-ts       same check for UPSTREAM_FRAC, using a reach that
                      splits flow between two downstream reaches
  sfr-man-pgap       a single-row PERIOD block for reach 2 only must
                      not misapply its value to reach 1
  sfr-ustrf-pgap     same check for UPSTREAM_FRAC
"""

import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "sfr-bedk",
    "sfr-bedkts",
    "sfr-rhkts",
    "sfr-rhkts-nstp",
    "sfr-manning-ts",
    "sfr-ustrf-ts",
    "sfr-man-pgap",
    "sfr-ustrf-pgap",
]


def _bedk_model(idx, ws, name):
    length_units = "meters"
    time_units = "seconds"

    nrow = 1
    ncol = 1
    nlay = 1
    delr = delc = 1.0

    nper = 4
    nstp = 2 if idx == 3 else 1
    tds_spd = [(1.0, nstp, 1.0)] * nper

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units=time_units,
        nper=nper,
        perioddata=tds_spd,
    )
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=-100.0,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,  # >0 means saturated thickness varies with computed head
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[((0, 0, 0), 1.0, 1e6)])

    # sfr data
    nreaches = 2

    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #        <man> <ncon> <ustrf> <ndv>
    if idx < 2:
        rhk1 = 0.0
        rhk2 = 0.0
    else:
        rhk1 = "bedk1"
        rhk2 = "bedk2"
    package_data = [
        (0, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, rhk1, 0.001, 0, 0.0, 0),
        (1, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, rhk2, 0.001, 0, 0.0, 0),
    ]
    connection_data = [
        (0,),
        (1,),
    ]

    if idx == 0:
        timeseries = False
        sfr_spd = {
            0: [
                (1, "bedk", 10.0),
            ],
            1: [
                (0, "bedk", 1.0),
                (1, "bedk", 5.0),
            ],
            2: [
                (0, "bedk", 5.0),
                (1, "bedk", 1.0),
            ],
            3: [
                (0, "bedk", 10.0),
                (1, "bedk", 0.0),
            ],
        }
    else:
        timeseries = True
        ts_names = ["bedk1", "bedk2"]
        ts_methods = ["linearend"] * len(ts_names)
        ts_data = [
            (0.0, 0.0, 10.0),
            (1.0, 0.0, 10.0),
            (2.0, 1.0, 5.0),
            (3.0, 5.0, 1.0),
            (4.0, 10.0, 0.0),
        ]
        if idx < 2:
            sfr_spd = {
                0: [
                    (0, "bedk", "bedk1"),
                    (1, "bedk", "bedk2"),
                ]
            }
        else:
            sfr_spd = None

    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwfr1", "sfr", (0,)),
            ("gwfr2", "sfr", (1,)),
        ],
        "filename": f"{name}.sfr.obs",
    }

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=1.0,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=sfr_spd,
        observations=sfr_obs,
        pname="SFR-1",
    )

    if timeseries:
        sfr.ts.initialize(
            filename=f"{name}.sfr.ts",
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=ts_methods,
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("head", "all"), ("budget", "all")],
    )

    return sim


def _manning_model(ws, name):
    # 2-reach chain; MANNING is PACKAGEDATA-TS-driven with no PERIOD
    # override, so any within-period change must come from per-timestep updates
    nper = 2
    nstp = 2
    tdis_rc = [(1.0, nstp, 1.0)] * nper

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-5, inner_dvclose=1e-6)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=0.0, botm=-100.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[((0, 0, 0), 1.0, 1e6)])

    packagedata = [
        (0, (0, 0, 0), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, "man1", 0, 0.0, 0),
        (1, (0, 0, 0), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, "man2", 0, 0.0, 0),
    ]
    connectiondata = [(0,), (1,)]
    ts_names = ["man1", "man2"]
    ts_data = [
        (0.0, 0.02, 0.02),
        (1.0, 0.02, 0.02),
        (2.0, 0.01, 0.03),
    ]
    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwfr1", "sfr", (0,)),
            ("gwfr2", "sfr", (1,)),
        ],
        "filename": f"{name}.sfr.obs",
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        length_conversion=1.0,
        time_conversion=1.0,
        nreaches=2,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [(0, "inflow", 1.0)]},
        observations=sfr_obs,
        pname="SFR-1",
    )
    sfr.ts.initialize(
        filename=f"{name}.sfr.ts",
        timeseries=ts_data,
        time_series_namerecord=ts_names,
        interpolation_methodrecord=["linearend"] * len(ts_names),
    )
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("head", "all"), ("budget", "all")])
    return sim


def _ustrf_model(ws, name):
    # 3-reach branch (reach 0 splits into 1/2); TS-driven ustrf1/ustrf2 must
    # stay complementary (sum to 1) every timestep, not just at period start
    nper = 2
    nstp = 2
    tdis_rc = [(1.0, nstp, 1.0)] * nper

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-5, inner_dvclose=1e-6)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=1.0, delc=1.0, top=0.0, botm=-100.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=[
            ((0, 0, 0), 1.0, 1e6),
            ((0, 0, 1), 1.0, 1e6),
            ((0, 0, 2), 1.0, 1e6),
        ],
    )

    packagedata = [
        (0, (0, 0, 0), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, 0.03, 2, 1.0, 0),
        (1, (0, 0, 1), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, 0.03, 1, "ustrf1", 0),
        (2, (0, 0, 2), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, 0.03, 1, "ustrf2", 0),
    ]
    connectiondata = [(0, -1, -2), (1, 0), (2, 0)]
    ts_data = [
        (0.0, 0.5, 0.5),
        (1.0, 0.5, 0.5),
        (2.0, 0.2, 0.8),
    ]
    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("usflow1", "upstream-flow", (1,)),
            ("usflow2", "upstream-flow", (2,)),
        ],
        "filename": f"{name}.sfr.obs",
    }
    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        length_conversion=1.0,
        time_conversion=1.0,
        nreaches=3,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [(0, "inflow", 1.0)]},
        observations=sfr_obs,
        pname="SFR-1",
    )
    sfr.ts.initialize(
        filename=f"{name}.sfr.ts",
        timeseries=ts_data,
        time_series_namerecord=["ustrf1", "ustrf2"],
        interpolation_methodrecord=["linearend", "linearend"],
    )
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("head", "all"), ("budget", "all")])
    return sim


def _get_periodgap_model(setting, value, ws, name):
    """A PERIOD row for reach 2 only must not misapply its value to
    reach 1's slot."""
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, outer_dvclose=1e-5, inner_dvclose=1e-6)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=0.0, botm=-100.0
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[((0, 0, 0), 1.0, 1e6)])

    packagedata = [
        (0, (0, 0, 0), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, 0.03, 0, 0.5, 0),
        (1, (0, 0, 0), 1.0, 1.0, 1e-3, 0.0, 1.0, 1e-5, 0.03, 0, 0.5, 0),
    ]
    connectiondata = [(0,), (1,)]
    flopy.mf6.ModflowGwfsfr(
        gwf,
        print_input=True,
        nreaches=2,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: [(1, setting, value)]},
        pname="SFR-1",
    )
    flopy.mf6.ModflowGwfoc(gwf, printrecord=[("budget", "all")])
    return sim


def _check_periodgap(keyword, value, test):
    """Reach 2's echoed value must match the PERIOD input, not reach 1's
    untouched default."""
    lst_path = test.workspace / f"{test.name}.lst"
    with open(lst_path) as f:
        text = f.read()
    match = re.search(rf"\n\s*2\s+{keyword}\s+([-\d.Ee+]+)", text)
    assert match, f"no {keyword} row found for reach 2 in the listing file"
    got = float(match.group(1))
    assert np.isclose(got, value), f"{keyword} for reach 2 expected {value}, got {got}"


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]
    if idx <= 3:
        return _bedk_model(idx, ws, name), None
    elif idx == 4:
        return _manning_model(ws, name), None
    elif idx == 5:
        return _ustrf_model(ws, name), None
    elif idx == 6:
        return _get_periodgap_model("manning", 0.5, ws, name), None
    else:
        return _get_periodgap_model("upstream_fraction", 0.9, ws, name), None


def check_output(idx, test):
    if idx == 6:
        _check_periodgap("MANNING", 0.5, test)
        return
    if idx == 7:
        _check_periodgap("UPSTREAM_FRAC", 0.9, test)
        return

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwf = sim.get_model()
    sfr = gwf.get_package("SFR-1")
    obs_data = sfr.output.obs().get_data()

    if idx <= 3:
        o1 = obs_data["GWFR1"]
        o2 = obs_data["GWFR2"][::-1]
        if idx == 3:
            # 4 periods x nstp=2 -> 8 records; bedk1 ramps in periods 3-4,
            # so GWFR1 must differ within each if it updates every timestep
            assert len(o1) == 8, (
                f"expected 8 obs records (4 periods x nstp=2), got {len(o1)}"
            )
            assert not np.isclose(o1[4], o1[5]), (
                f"GWFR1 unchanged in period 3 (t=2.5->3.0): "
                f"{o1[4]:.6g} == {o1[5]:.6g}; rhk_src re-sync in sfr_ad "
                "may have failed"
            )
            assert not np.isclose(o1[6], o1[7]), (
                f"GWFR1 unchanged in period 4 (t=3.5->4.0): "
                f"{o1[6]:.6g} == {o1[7]:.6g}; rhk_src re-sync in sfr_ad "
                "may have failed"
            )
        else:
            assert np.allclose(o1, o2), (
                f"GWFR1 ({o1}) not equal to reversed GWFR2 ({o2})"
            )
    elif idx == 4:
        gwfr1 = obs_data["GWFR1"]
        # 2 periods x nstp=2 -> 4 records; MANNING ramps in period 1, so
        # GWFR1 must differ within it if the value updates every timestep
        assert len(gwfr1) == 4, f"expected 4 obs records, got {len(gwfr1)}"
        assert not np.isclose(gwfr1[2], gwfr1[3]), (
            f"GWFR1 unchanged within period 1 (t=1.5->2.0): "
            f"{gwfr1[2]:.6g} == {gwfr1[3]:.6g}; MANNING_SET baseline "
            "resync in sfr_ad may have failed"
        )
    elif idx == 5:
        usflow1 = obs_data["USFLOW1"]
        # 2 periods x nstp=2 -> 4 records; USTRF1 ramps in period 1, so
        # its upstream-flow share must differ within it, not just once
        assert len(usflow1) == 4, f"expected 4 obs records, got {len(usflow1)}"
        assert not np.isclose(usflow1[2], usflow1[3]), (
            f"USFLOW1 unchanged within period 1 (t=1.5->2.0): "
            f"{usflow1[2]:.6g} == {usflow1[3]:.6g}; USTRF_SET baseline "
            "resync in sfr_ad may have failed"
        )


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

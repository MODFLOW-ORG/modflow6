"""LAK lakes with no GWF connections (NLAKECONN = 0).

A lake with no connections to the aquifer (e.g. a reservoir) is supported
if it has a lake table, which defines its bottom, top, volume and surface
area. Its stage is determined by its own budget: inflows, outlets,
rainfall, evaporation, withdrawal and storage.

- lak-noconn: lake 1 is connected to the aquifer and drains through a
  SPECIFIED outlet into lake 2, which has no connections. Lake 2 also
  receives rainfall over the table surface area.
- lak-noconn-imp: as above, with the IMPLICIT formulation.
- lak-noconn-only: the only lake has no connections, so there is no
  CONNECTIONDATA block; it is filled by a specified inflow and rainfall.

Lake 2's table has a constant surface area, so its volume is linear in
stage and its stage at the end of each time step can be computed exactly
from the cumulative inflow.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-noconn", "lak-noconn-imp", "lak-noconn-only"]

nper = 1
nstp = 5
perlen = 5.0
area = 1000.0  # table surface area of the unconnected lake
strt = 1.0  # starting stage of the unconnected lake
qin = 100.0  # outlet or external inflow to the unconnected lake
rain = 0.01  # rainfall rate on the unconnected lake


def build_models(idx, test):
    name = cases[idx]
    implicit = "imp" in name
    only = "only" in name

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, nstp, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, sy=0.1, ss=1e-5, iconvert=1, transient={0: True})
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(0, 0, 0), 5.0], [(0, 0, 2), 5.0]]
    )

    tab_fname = f"{name}.laktab"
    if only:
        # one lake, no connections, no CONNECTIONDATA block
        ilak = 0
        packagedata = [(0, strt, 0)]
        connectiondata = None
        noutlets = 0
        outlets = None
        perioddata = [
            (0, "status", "active"),
            (0, "inflow", qin),
            (0, "rainfall", rain),
        ]
    else:
        # lake 1 (connected) drains to lake 2 (unconnected)
        ilak = 1
        packagedata = [(0, 11.0, 1), (1, strt, 0)]
        connectiondata = [(0, 0, (0, 0, 1), "vertical", 0.1, 0.0, 0.0, 0.0, 0.0)]
        noutlets = 1
        outlets = [(0, 0, 1, "specified", 0.0, 0.0, 0.0, 0.0)]
        perioddata = [
            (0, "status", "active"),
            (1, "status", "active"),
            (0, "inflow", 500.0),
            (0, "rate", -qin),
            (1, "rainfall", rain),
        ]

    lak = flopy.mf6.ModflowGwflak(
        gwf,
        implicit=implicit,
        print_stage=True,
        print_flows=True,
        save_flows=True,
        nlakes=len(packagedata),
        noutlets=noutlets,
        ntables=1,
        tables=[(ilak, tab_fname)],
        packagedata=packagedata,
        connectiondata=connectiondata,
        outlets=outlets,
        perioddata=perioddata,
        budget_filerecord=f"{name}.lak.cbc",
        pname="lak-1",
    )
    stages = np.linspace(0.0, 10.0, 11)
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(stages),
        ncol=3,
        table=[(s, area * s, area) for s in stages],
        filename=tab_fname,
        parent_file=lak,
    )
    lak.obs.initialize(
        filename=f"{name}.lak.obs",
        digits=15,
        continuous={
            f"{name}.lak.obs.csv": [
                ("stage", "stage", (ilak,)),
                ("volume", "volume", (ilak,)),
                ("rain", "rainfall", (ilak,)),
                ("storage", "storage", (ilak,)),
            ]
        },
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def check_output(idx, test):
    name = cases[idx]
    obs = np.genfromtxt(
        test.workspace / f"{name}.lak.obs.csv", delimiter=",", names=True
    )
    dt = perlen / nstp
    t = dt * np.arange(1, nstp + 1)
    expected_vol = area * strt + (qin + rain * area) * t
    expected_stage = expected_vol / area
    assert np.allclose(obs["time"], t)
    assert np.allclose(obs["VOLUME"], expected_vol, atol=1e-6), (
        obs["VOLUME"],
        expected_vol,
    )
    assert np.allclose(obs["STAGE"], expected_stage, atol=1e-9), (
        obs["STAGE"],
        expected_stage,
    )
    # rainfall is applied over the table surface area
    assert np.allclose(obs["RAIN"], rain * area)
    # storage takes up all of the inflow (negative = into storage)
    assert np.allclose(obs["STORAGE"], -(qin + rain * area), atol=1e-6)


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

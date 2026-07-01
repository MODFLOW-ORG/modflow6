"""
Test that stress packages issue a warning (not an error) when a boundary
condition entry references a cell that is not in the active grid domain.
The model should run to completion and the listing file should contain a
warning about the inactive cell.

Packages tested: RCH, WEL, GHB, DRN, RIV, CHD
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "inact_rch",
    "inact_wel",
    "inact_ghb",
    "inact_drn",
    "inact_riv",
    "inact_chd",
]

nlay, nrow, ncol = 1, 3, 3
inactive_cell = (0, 0, 0)
chd_cell = (0, 2, 2)
top = 10.0
botm = 0.0
strt = 5.0
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[inactive_cell] = 0


def _base_sim(test, name):
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[chd_cell, strt]]},
    )

    return sim, gwf


def build_models(idx, test):
    name = cases[idx]
    sim, gwf = _base_sim(test, name)

    if idx == 0:
        # RCH: one entry on the inactive cell, one on an active cell
        flopy.mf6.ModflowGwfrch(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, 1e-3],
                    [(0, 1, 1), 1e-3],
                ]
            },
        )

    elif idx == 1:
        # WEL: one entry on the inactive cell
        flopy.mf6.ModflowGwfwel(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, -1.0],
                    [(0, 1, 1), -1.0],
                ]
            },
        )

    elif idx == 2:
        # GHB: one entry on the inactive cell
        flopy.mf6.ModflowGwfghb(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, strt, 1.0],
                    [(0, 1, 1), strt, 1.0],
                ]
            },
        )

    elif idx == 3:
        # DRN: one entry on the inactive cell
        flopy.mf6.ModflowGwfdrn(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, 3.0, 1.0],
                    [(0, 1, 1), 3.0, 1.0],
                ]
            },
        )

    elif idx == 4:
        # RIV: one entry on the inactive cell
        flopy.mf6.ModflowGwfriv(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, strt, 1.0, botm],
                    [(0, 1, 1), strt, 1.0, botm],
                ]
            },
        )

    elif idx == 5:
        # CHD: one entry on the inactive cell (the base sim already has a
        # CHD on an active cell; add a separate package with the bad entry)
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={
                0: [
                    [inactive_cell, strt],
                    [(0, 0, 1), strt],
                ]
            },
            pname="chd_bad",
            filename=f"{name}_bad.chd",
        )

    return sim, None


def check_output(idx, test):
    mfsim_lst = os.path.join(test.workspace, "mfsim.lst")
    assert os.path.isfile(mfsim_lst)
    lst_text = open(mfsim_lst).read()
    assert "Normal termination" in lst_text
    warning_phrase = "outside active grid domain"
    assert warning_phrase.lower() in lst_text.lower()


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

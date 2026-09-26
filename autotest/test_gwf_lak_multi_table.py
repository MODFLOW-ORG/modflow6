"""Two lakes with different TAB6 tables (differing NROW/NCOL, one with
BAREA) listed out of IFNO order, verifying per-lake table matching by
IFNO rather than TABLES-block position. Since both lakes are STATUS
CONSTANT, expected VOLUME/SURFACE-AREA are computed by interpolating
each lake's own table at its own STRT.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-multi-table"]

nlay, nrow, ncol = 1, 1, 2
top = 10.0
botm = -10.0
k11 = 10.0

# lake 0: HORIZONTAL connection, 3-column table (no BAREA), 3 rows
lak0_strt = 7.0
lak0_bedleak = 1.0
lak0_table = [
    [0.0, 0.0, 10.0],
    [5.0, 50.0, 10.0],
    [10.0, 150.0, 15.0],
]

# lake 1: EMBEDDEDH connection, 4-column table (BAREA present), 5 rows --
# deliberately a different row count and column count than lake 0's table
lak1_strt = 5.0
lak1_bedleak = 1.0
lak1_table = [
    [0.0, 0.0, 10.0, 10.0],
    [2.0, 20.0, 10.0, 10.0],
    [4.0, 45.0, 12.0, 11.0],
    [6.0, 75.0, 14.0, 12.0],
    [8.0, 120.0, 16.0, 13.0],
]


def _interp(table, stage, col):
    """Linear interpolation matching MF6's own lak_linear_interpolation:
    table rows are [stage, volume, sarea, (barea)], col selects the
    output column (1=volume, 2=sarea)."""
    stages = [row[0] for row in table]
    for i in range(len(stages) - 1):
        if stages[i] <= stage <= stages[i + 1]:
            frac = (stage - stages[i]) / (stages[i + 1] - stages[i])
            v0 = table[i][col]
            v1 = table[i + 1][col]
            return v0 + frac * (v1 - v0)
    raise ValueError(f"stage {stage} outside table range {stages}")


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, print_option="SUMMARY", complexity="simple")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k11)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)

    lak_packagedata = [
        [0, lak0_strt, 1],
        [1, lak1_strt, 1],
    ]
    lak_conn = [
        [0, 0, (0, 0, 0), "horizontal", lak0_bedleak, -5.0, 5.0, 10.0, 10.0],
        [1, 0, (0, 0, 1), "embeddedh", lak1_bedleak, 0.0, 0.0, 1.0, 0.0],
    ]
    lak_spd = {
        0: [
            [0, "status", "constant"],
            [1, "status", "constant"],
        ]
    }

    tab0_filename = f"{name}.lak0.tab"
    tab1_filename = f"{name}.lak1.tab"
    budpth = f"{name}.lak.cbc"
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        save_flows=True,
        print_stage=True,
        nlakes=2,
        noutlets=0,
        ntables=2,
        tables=[[1, tab1_filename], [0, tab0_filename]],  # out-of-order IFNO
        packagedata=lak_packagedata,
        connectiondata=lak_conn,
        perioddata=lak_spd,
        budget_filerecord=budpth,
        pname="LAK-1",
    )
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(lak0_table),
        ncol=len(lak0_table[0]),
        table=lak0_table,
        filename=tab0_filename,
        pname="LAK0TAB",
        parent_file=lak,
    )
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(lak1_table),
        ncol=len(lak1_table[0]),
        table=lak1_table,
        filename=tab1_filename,
        pname="LAK1TAB",
        parent_file=lak,
    )

    obs_file = f"{name}.lak.obs"
    csv_file = obs_file + ".csv"
    obs_dict = {
        csv_file: [
            ("stage0", "stage", (0,)),
            ("stage1", "stage", (1,)),
            ("volume0", "volume", (0,)),
            ("volume1", "volume", (1,)),
            ("sarea0", "surface-area", (0,)),
            ("sarea1", "surface-area", (1,)),
        ]
    }
    lak.obs.initialize(
        filename=obs_file, digits=10, print_input=True, continuous=obs_dict
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
    csv_file = os.path.join(test.workspace, f"{name}.lak.obs.csv")
    obs = np.atleast_1d(np.genfromtxt(csv_file, names=True, delimiter=","))

    stage0 = obs["STAGE0"][-1]
    stage1 = obs["STAGE1"][-1]
    assert np.isclose(stage0, lak0_strt), f"lake 0 stage {stage0} != STRT {lak0_strt}"
    assert np.isclose(stage1, lak1_strt), f"lake 1 stage {stage1} != STRT {lak1_strt}"

    expected_vol0 = _interp(lak0_table, lak0_strt, 1)
    expected_sarea0 = _interp(lak0_table, lak0_strt, 2)
    expected_vol1 = _interp(lak1_table, lak1_strt, 1)
    expected_sarea1 = _interp(lak1_table, lak1_strt, 2)

    vol0 = obs["VOLUME0"][-1]
    sarea0 = obs["SAREA0"][-1]
    vol1 = obs["VOLUME1"][-1]
    sarea1 = obs["SAREA1"][-1]

    assert np.isclose(vol0, expected_vol0), (
        f"lake 0 volume {vol0} != expected {expected_vol0} "
        "(wrong table applied to lake 0?)"
    )
    assert np.isclose(sarea0, expected_sarea0), (
        f"lake 0 surface-area {sarea0} != expected {expected_sarea0}"
    )
    assert np.isclose(vol1, expected_vol1), (
        f"lake 1 volume {vol1} != expected {expected_vol1} "
        "(wrong table applied to lake 1, or row-count offset bug?)"
    )
    assert np.isclose(sarea1, expected_sarea1), (
        f"lake 1 surface-area {sarea1} != expected {expected_sarea1}"
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

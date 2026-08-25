"""
These test cases exercise cell "snap" or "clamp" behavior in lateral and
vertical dimensions. There is currently no distance limit; particles may
begin arbitrarily far outside the cell.

For rectilinear cells using Pollock's method, a particle will be clamped
to the point within the cell nearest to the point outside: a "box clamp".

For cells using the generalized (ternary) tracking method, a particle is
not guaranteed to end up at the nearest point within the cell.
"""

import flopy
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

DISTANCES = {
    # comparable to the padding clamp_bary() applies (DSAME * DEP3)
    "t": 1e-7,
    # a fraction of a cell width
    "m": 0.4,
    # several cells away, similar in scale to the mismatch case above
    "f": 4.5,
}
METHODS = {"p": "pollock", "n": "ternary"}
DIRECTIONS = {"": "x", "z": "z"}

cases = [
    (f"{dircode}{mcode}{dcode}", method, offset, direction)
    for dircode, direction in DIRECTIONS.items()
    for mcode, method in METHODS.items()
    for dcode, offset in DISTANCES.items()
]


def get_gridprops_rect():
    return get_disv_kwargs(
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
        delr=1.0,
        delc=1.0,
        tp=FlopyReadmeCase.top,
        botm=FlopyReadmeCase.botm,
    )


def build_gwf_sim_disv(name, ws, mf6):
    gridprops = get_gridprops_rect()

    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name=mf6, version="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

    gwfname = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    flopy.mf6.ModflowGwfdisv(gwf, **gridprops)
    flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, save_saturation=True)
    flopy.mf6.ModflowGwfic(gwf, strt=FlopyReadmeCase.top)

    ncpl = gridprops["ncpl"]
    chd_data = [[(0, 0), FlopyReadmeCase.top], [(0, ncpl - 1), 0.0]]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_data)

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowIms(sim)
    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6, method, offset, direction):
    gridprops = get_gridprops_rect()

    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, version="mf6", sim_ws=prt_ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowPrtdisv(prt, **gridprops)
    flopy.mf6.ModflowPrtmip(prt, porosity=FlopyReadmeCase.porosity)

    if direction == "x":
        release_x = 1.0 + offset
        release_y = 9.5
        release_z = 0.5
    else:
        release_x = 0.5
        release_y = 9.5
        release_z = FlopyReadmeCase.botm[0] - offset
    releasepts = [[0, (0, 0), release_x, release_y, release_z]]

    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: [("FIRST",)]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        coordinate_check_method="none",
        dev_forceternary=method == "ternary",
        print_input=True,
        extend_tracking=True,
    )

    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    gwf_name = get_model_name(name, "gwf")
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(test, method, offset, direction):
    gwf_sim = build_gwf_sim_disv(test.name, test.workspace, test.targets["mf6"])
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets["mf6"],
        method,
        offset,
        direction,
    )
    return gwf_sim, prt_sim


def check_output(test, offset, direction):
    prt_name = get_model_name(test.name, "prt")
    prt_ws = test.workspace / "prt"
    mf6_pls = pd.read_csv(prt_ws / f"{prt_name}.trk.csv", na_filter=False)

    assert len(mf6_pls) > 1

    # the release event should report the raw release coordinates
    release = mf6_pls.iloc[0]
    if direction == "x":
        assert release["x"] == pytest.approx(1.0 + offset)
        assert release["y"] == pytest.approx(9.5)
    else:
        assert release["x"] == pytest.approx(0.5)
        assert release["y"] == pytest.approx(9.5)
        assert release["z"] == pytest.approx(FlopyReadmeCase.botm[0] - offset)

    # subsequent events should report the clamped coordinates
    ncol = FlopyReadmeCase.ncol
    nrow = FlopyReadmeCase.nrow
    top = FlopyReadmeCase.top
    tracking = mf6_pls.iloc[1:]
    assert tracking["x"].between(-1e-6, ncol + 1e-6).all()
    assert tracking["y"].between(-1e-6, nrow + 1e-6).all()
    assert tracking["z"].between(-1e-6, top + 1e-6).all()

    # time should never decrease (texit should never be negative)
    assert (mf6_pls["t"].diff().dropna() >= -1e-9).all()


@pytest.mark.developmode
@pytest.mark.parametrize("name, method, offset, direction", cases)
def test_mf6model(name, method, offset, direction, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, method, offset, direction),
        check=lambda t: check_output(t, offset, direction),
        targets=targets,
        compare=None,
    )
    test.run()

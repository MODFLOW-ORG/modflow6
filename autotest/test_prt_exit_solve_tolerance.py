"""
Test exit_solve_tolerance validation and warning behavior.

Tests four scenarios:
1. High tolerance (> 0.01): Should issue a warning showing the actual value
2. Normal tolerance (<= 0.01): Should run without warning
3. Negative tolerance (< 0): Should error and fail
4. Zero tolerance (= 0): Should error and fail
"""

import os

import flopy
import pytest
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtex"
cases = [
    f"{simname}_high",  # 0.1, should warn
    f"{simname}_ok",  # 0.005, should not warn
    f"{simname}_neg",  # -0.1, should error
    f"{simname}_zero",  # 0.0, should error
]


def get_tolerance(name):
    """Get the exit_solve_tolerance value for the test case."""
    if "high" in name:
        return 0.1
    elif "ok" in name:
        return 0.005
    elif "neg" in name:
        return -0.1
    elif "zero" in name:
        return 0.0
    raise ValueError(f"Unknown test case: {name}")


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # create tdis package
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=FlopyReadmeCase.nper,
        perioddata=[
            (
                FlopyReadmeCase.perlen,
                FlopyReadmeCase.nstp,
                FlopyReadmeCase.tsmult,
            )
        ],
    )

    # create prt model
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    # create prt discretization
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # create a single release point
    releasepts = [[0, 0, 0, 0, 0.5, 0.5, 0.5]]

    # create prp package with the specified exit_solve_tolerance
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"

    tolerance = get_tolerance(name)

    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: [("FIRST",)]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        exit_solve_tolerance=tolerance,
        coordinate_check_method="none",  # Disable coordinate checking for this test
        print_input=True,
        extend_tracking=True,
    )

    # create output control package
    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    # create the flow model interface
    gwf_name = get_model_name(name, "gwf")
    rel_gwf_folder = os.path.relpath(gwf_ws, start=prt_ws)
    gwf_budget_file = f"{rel_gwf_folder}/{gwf_name}.bud"
    gwf_head_file = f"{rel_gwf_folder}/{gwf_name}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    # add explicit model solution
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(test):
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace / "gwf", test.targets["mf6"]
    )
    prt_sim = build_prt_sim(
        test.name,
        test.workspace / "gwf",
        test.workspace / "prt",
        test.targets["mf6"],
    )
    return gwf_sim, prt_sim


def check_output(test):
    name = test.name
    tolerance = get_tolerance(name)

    # For the PRT simulation (second simulation), check output buffer
    # test.buffs[0] is GWF, test.buffs[1] is PRT
    buff = test.buffs[1]

    if "high" in name:
        # Should have warning message for high tolerance (> 0.01)
        assert any("WARNING: EXIT_SOLVE_TOLERANCE is set to" in l for l in buff), (
            "Expected warning for high exit_solve_tolerance not found"
        )

        # Check that the warning shows the actual tolerance value
        buff_text = " ".join(buff)
        assert "0.100" in buff_text or "0.1" in buff_text or "1.00E-01" in buff_text, (
            "Warning should show tolerance value 0.1 in the output"
        )

        # Check that it mentions the default value
        assert "default value" in buff_text, "Warning should mention the default value"

    elif "ok" in name:
        # Should NOT have warning message for tolerance <= 0.01
        assert not any("WARNING: EXIT_SOLVE_TOLERANCE is set to" in l for l in buff), (
            "Unexpected warning for acceptable exit_solve_tolerance"
        )

    elif "neg" in name or "zero" in name:
        # Should have error message for negative or zero tolerance
        assert any("EXIT_SOLVE_TOLERANCE MUST BE POSITIVE" in l for l in buff), (
            "Expected error for non-positive exit_solve_tolerance not found"
        )


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build_models,
        check=check_output,
        targets=targets,
        compare=None,
        # For non-positive tolerance cases: GWF should succeed, PRT should fail
        xfail=[False, True] if ("neg" in name or "zero" in name) else False,
    )
    test.run()

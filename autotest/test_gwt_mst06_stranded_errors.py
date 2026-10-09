"""
Tests the errors and warnings of the MST stranded mass option, using the single
cell of test_gwt_mst06_stranded_cell.py.

Cases:
  - ist               : STRANDED_MASS with an active IST Package is an error
                        (xfail).
  - fileout           : STRANDED FILEOUT without STRANDED_MASS is an error
                        (xfail).
  - zero_order        : zero-order decay runs and warns that stranded mass does
                        not decay.
  - sy_above_porosity : a specific yield greater than the porosity runs and warns
                        that no solute is held back from residual water.
  - fmi_first_step    : a transport model that reads its flows from files, with
                        a water table that moves in the first time step, runs,
                        warns that no mass was stranded then, and strands none;
                        so does a cell that the standard formulation dries in
                        the first time step, for which the flow model saves no
                        storage; a leading period in which the water table is
                        still, or a cell whose head stays above its top while
                        specific storage releases water, writes no warning.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework
from test_gwt_mst06_stranded_cell import cinit, get_model, paper, porosity


def listing_text(ws):
    """Simulation listing with line breaks removed, for wrapped messages."""
    return " ".join((ws / "mfsim.lst").read_text().split())


def run_framework(function_tmpdir, targets, build, check, xfail=False):
    TestFramework(
        name="cell",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        xfail=xfail,
    ).run()


def test_ist(function_tmpdir, targets):
    def build(test):
        sim = get_model(test.workspace, paper, 1, porosity, True)
        gwt = sim.get_model("gwt")
        flopy.mf6.ModflowGwtist(gwt, porosity=0.05, volfrac=0.1, zetaim=0.1)
        return sim

    def check(test):
        assert (
            "STRANDED_MASS is specified in the MST package, but an immobile "
            "domain (IST) is active" in listing_text(test.workspace)
        )

    run_framework(function_tmpdir, targets, build, check, xfail=True)


def test_fileout(function_tmpdir, targets):
    def build(test):
        return get_model(
            test.workspace,
            paper,
            1,
            porosity,
            True,
            strand=False,
            mst_kwargs={"stranded_filerecord": [("gwt.strand.bin",)]},
        )

    def check(test):
        assert (
            "STRANDED FILEOUT was specified but the STRANDED_MASS keyword was not "
            "specified" in listing_text(test.workspace)
        )

    run_framework(function_tmpdir, targets, build, check, xfail=True)


def test_zero_order(function_tmpdir, targets):
    def build(test):
        return get_model(
            test.workspace,
            paper,
            1,
            porosity,
            True,
            mst_kwargs={
                "zero_order_decay": True,
                "decay": 1.0e-3,
                "decay_sorbed": 1.0e-6,
            },
        )

    def check(test):
        assert "Zero-order decay is active with STRANDED_MASS" in listing_text(
            test.workspace
        )

    run_framework(function_tmpdir, targets, build, check)


def test_sy_above_porosity(function_tmpdir, targets):
    def build(test):
        return get_model(test.workspace, paper, 1, 1.5 * porosity, True)

    def check(test):
        assert "Specific yield is greater than the mobile porosity" in listing_text(
            test.workspace
        )

    run_framework(function_tmpdir, targets, build, check)


# the first stress period of each case, the starting head, the specific
# storage, and whether the flow model uses the Newton formulation; only in the
# moved and dry cases does the water table move in the first step
first_step_cases = {
    "moved": ([], 10.0, 0.0, True),
    # the standard formulation dries the cell in the first time step
    "dry": ([], 10.0, 0.0, False),
    "lead": ([(1.0, [(0.0, 0.0)])], 10.0, 0.0, True),
    # 1 m3/d lowers the head from 15 to 14 m through specific storage, above
    # the 10 m top of the cell
    "confined": ([(1.0, [(-1.0, 0.0)])], 15.0, 1.0e-3, True),
}


@pytest.mark.parametrize("case", list(first_step_cases))
def test_fmi_first_step(function_tmpdir, targets, case):
    """A warning is written only when the water table moves in the first step."""
    lead, strt, ss, newton = first_step_cases[case]
    exe = str(targets["mf6"])
    periods = lead + paper

    # the flow model alone, writing the head and budget files
    flow_ws = function_tmpdir / "flow"
    sim = get_model(
        flow_ws,
        periods,
        1,
        porosity,
        True,
        strt=strt,
        ss=ss,
        transport=False,
        newton=newton,
    )
    sim.exe_name = exe
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"flow simulation failed\n{buff}"

    # the transport model alone, reading the flows through FMI
    ws = function_tmpdir / "transport"
    sim = flopy.mf6.MFSimulation(sim_name="t", sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=len(periods),
        perioddata=[(perlen, 1, 1.0) for perlen, _ in periods],
    )
    flopy.mf6.ModflowIms(
        sim, complexity="moderate", linear_acceleration="bicgstab", filename="t.ims"
    )
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    flopy.mf6.ModflowGwtdis(gwt, nrow=1, ncol=1, delr=10.0, delc=10.0, top=10.0)
    flopy.mf6.ModflowGwtic(gwt, strt=cinit)
    flopy.mf6.ModflowGwtmst(
        gwt,
        porosity=porosity,
        sorption="linear",
        bulk_density=1500.0,
        distcoef=1.0e-4,
        stranded_mass=True,
        stranded_filerecord=[("gwt.strand.bin",)],
    )
    flopy.mf6.ModflowGwtfmi(
        gwt,
        packagedata=[
            ("GWFHEAD", str(flow_ws / "gwf.hds"), None),
            ("GWFBUDGET", str(flow_ws / "gwf.cbc"), None),
        ],
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    # stranded mass is written with the concentration
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"transport simulation failed\n{buff}"

    text = listing_text(ws)
    warned = (
        "The water table moved during the first time step" in text
        or "At least one cell was dry at the end of the first time step" in text
    )
    assert warned == (case in ("moved", "dry")), (
        f"a warning was expected only when the water table moved, {case=} {warned=}"
    )
    # the saturation the first time step starts from is not in the budget
    # file, so no mass is stranded during it; the leading period makes the
    # drainage happen in the second time step, where mass is stranded
    sf = flopy.utils.HeadFile(ws / "gwt.strand.bin", text="STRANDED")
    stranded = np.array([sf.get_data(totim=t).ravel()[0] for t in sf.get_times()])
    assert stranded[0] == 0.0, f"mass was stranded in the first time step: {stranded}"
    if case == "lead":
        assert stranded[1] > 0.0, f"no mass was stranded after the lead: {stranded}"

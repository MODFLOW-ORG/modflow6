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
                        a water table that moves in the first time step, runs and
                        warns that no mass was stranded then.
"""

import flopy
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


@pytest.mark.parametrize("lead", [False, True])
def test_fmi_first_step(function_tmpdir, targets, lead):
    """A leading stress period in which the water table is still avoids it."""
    exe = str(targets["mf6"])
    periods = [(1.0, [(0.0, 0.0)])] + paper if lead else paper

    # the flow model alone, writing the head and budget files
    flow_ws = function_tmpdir / "flow"
    sim = get_model(flow_ws, periods, 1, porosity, True, transport=False)
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
    )
    flopy.mf6.ModflowGwtfmi(
        gwt,
        packagedata=[
            ("GWFHEAD", str(flow_ws / "gwf.hds"), None),
            ("GWFBUDGET", str(flow_ws / "gwf.cbc"), None),
        ],
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[["WEL-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(gwt, concentration_filerecord="gwt.ucn")
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True)
    assert success, f"transport simulation failed\n{buff}"

    warned = "The water table moved during the first time step" in listing_text(ws)
    assert warned != lead, (
        "a warning was expected only without the leading stress period, "
        f"lead={lead} warned={warned}"
    )

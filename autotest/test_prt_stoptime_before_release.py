"""
Reproduce https://github.com/MODFLOW-ORG/modflow6/issues/2941

If a PRP package's STOPTIME precedes a particle's release time, the particle
should not be released, but terminate immediately (status 8), with a warning.

Two configurations are tested, matching the two models attached to the issue:
    - rt: an explicit RELEASETIMES entry falls after STOPTIME
    - pd: a PERIOD block release (FIRST) fills forward into a later stress
      period whose start time falls after STOPTIME
"""

import flopy
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtstpb4"
cases = [f"{simname}rt", f"{simname}pd"]

# release point matching cell (0, 0, 0) in the FlopyReadmeCase grid
releasepts = [[0, 0, 0, 0, 0.1, 9.1, 0.5]]


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    # create simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    # two stress periods so a later-period release time/step can fall
    # after STOPTIME, which is set to expire during the first period
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=2,
        perioddata=[
            (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult),
            (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult),
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
        top=FlopyReadmeCase.top,
        botm=FlopyReadmeCase.botm,
    )

    # create mip package
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)

    # create prp package
    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"

    if name.endswith("rt"):
        # release time (1.5) falls in period 2, after stoptime (0.5)
        perioddata = None
        nreleasetimes = 1
        releasetimes = [(1.5,)]
    else:
        # FIRST in period 1 releases at t=0.0 (before stoptime), but fills
        # forward into period 2, releasing again at t=1.0 (after stoptime)
        perioddata = {0: [("FIRST",)]}
        nreleasetimes = None
        releasetimes = None

    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata=perioddata,
        nreleasetimes=nreleasetimes,
        releasetimes=releasetimes,
        stoptime=0.5,
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
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
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
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
        test.name, test.workspace, test.targets["mf6"]
    )
    # GWF sim also needs 2 stress periods to match the PRT model's TDIS
    tdis = gwf_sim.get_package("tdis")
    tdis.nper = 2
    tdis.perioddata = [
        (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult),
        (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult),
    ]
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets["mf6"],
    )
    return gwf_sim, prt_sim


def check_output(test):
    # expect warning that particles are unreleased as stop time < release time
    lst = " ".join((test.workspace / "prt" / "mfsim.lst").read_text().split())
    assert "particle will not be released" in lst

    # particles scheduled after the package's stop time (0.5) should have just
    # a single event, termination with status 8 (permanently unreleased)
    prt_name = get_model_name(test.name, "prt")
    trk = pd.read_csv(test.workspace / "prt" / f"{prt_name}.prp.trk.csv")
    late = trk[trk["trelease"] > 0.5]
    assert len(late) == 1
    assert late.iloc[0]["istatus"] == 8


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=build_models,
        check=check_output,
        targets=targets,
        compare=None,
    )
    test.run()

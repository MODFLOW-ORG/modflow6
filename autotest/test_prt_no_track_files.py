"""
A PRT model without track output (no TRACK or TRACKCSV FILEOUT in OC
or any PRP) should run to normal termination with a warning. This
used to segfault on the first time step when flushing the event buffer.

Two cases: events buffered in memory (default) and in a scratch file.
"""

import flopy
import pytest
from framework import TestFramework
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtnotrk"
cases = [f"{simname}mem", f"{simname}scr"]


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )
    flopy.mf6.modflow.mftdis.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=FlopyReadmeCase.nper,
        perioddata=[
            (FlopyReadmeCase.perlen, FlopyReadmeCase.nstp, FlopyReadmeCase.tsmult)
        ],
    )

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
        top=FlopyReadmeCase.top,
        botm=FlopyReadmeCase.botm,
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=FlopyReadmeCase.porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(FlopyReadmeCase.releasepts_prt),
        packagedata=FlopyReadmeCase.releasepts_prt,
        perioddata={0: ["FIRST"]},
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=[f"{prt_name}.cbc"],
        saverecord=[("BUDGET", "ALL")],
        scratch_buffer=name.endswith("scr"),
    )

    gwf_name = get_model_name(name, "gwf")
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_ws / f"{gwf_name}.hds"),
            ("GWFBUDGET", gwf_ws / f"{gwf_name}.bud"),
        ],
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(test):
    gwf_sim = FlopyReadmeCase.get_gwf_sim(
        test.name, test.workspace, test.targets["mf6"]
    )
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets["mf6"],
    )
    return gwf_sim, prt_sim


def check_output(test):
    prt_ws = test.workspace / "prt"
    lst = " ".join((prt_ws / "mfsim.lst").read_text().split())
    assert "No particle track output configured" in lst
    assert not list(prt_ws.glob("*.trk*"))

    prt_name = get_model_name(test.name, "prt")
    assert (prt_ws / f"{prt_name}.cbc").is_file()


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

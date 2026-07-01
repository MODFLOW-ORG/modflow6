"""
Bigger PRT benchmark on the Voronoi grid from test_prt_voronoi1.py.

Uses the same simple left-to-right flow field but with many more particles,
generating many more track events. Meant to benchmark event buffer options,
in-memory vs scratch file, at a scale where I/O overhead differences start
to be measurable; very likely larger than configured here, but the runtime
also grows quickly, so tune particle count and TDIS granularity as needed.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from flopy.discretization import VertexGrid
from framework import TestFramework
from prt_test_utils import get_model_name
from test_prt_voronoi1 import botm, build_gwf_sim, get_grid, nlay, porosity, top

simname = "prtvorperf"
ntracktimes = 500
tracktimes = list(np.linspace(0, 40000, ntracktimes))
rpts = [[20, y, 0.5] for y in np.linspace(0.5, 999.5, 10000)]


def build_prt_sim(name, gwf_ws, prt_ws, targets, scratch_buffer=False):
    prt_ws = Path(prt_ws)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")

    grid = get_grid(prt_ws / "grid", targets)
    vgrid = VertexGrid(**grid.get_gridprops_vertexgrid(), nlay=1)

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=targets["mf6"], sim_ws=prt_ws
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", perioddata=[[1.0, 1, 1.0]])
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdisv(
        prt, nlay=nlay, **grid.get_disv_gridprops(), top=top, botm=botm
    )
    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    prpdata = [
        (i, (0, vgrid.intersect(p[0], p[1])), p[0], p[1], p[2])
        for i, p in enumerate(rpts)
    ]
    prp_track_file = f"{prt_name}.prp.trk"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        boundnames=True,
        stop_at_weak_sink=True,
        exit_solve_tolerance=1e-10,
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_exit=True,
        track_release=True,
        track_terminate=True,
        track_usertime=True,
        ntracktimes=ntracktimes,
        tracktimes=[(t,) for t in tracktimes],
        scratch_buffer=scratch_buffer,
    )
    gwf_budget_file = gwf_ws / f"{gwf_name}.bud"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    flopy.mf6.ModflowPrtfmi(
        prt, packagedata=[("GWFHEAD", gwf_head_file), ("GWFBUDGET", gwf_budget_file)]
    )
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(test, scratch_buffer=False):
    gwf_sim, _ = build_gwf_sim(test.name, test.workspace, test.targets)
    prt_sim = build_prt_sim(
        test.name,
        test.workspace,
        test.workspace / "prt",
        test.targets,
        scratch_buffer=scratch_buffer,
    )
    return gwf_sim, prt_sim


def check_output(test):
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(test.name, "prt")
    trk = prt_ws / f"{prt_name}.prp.trk"
    assert trk.exists() and trk.stat().st_size > 0


@pytest.mark.skip(reason="run manually")  # uncomment to run
@pytest.mark.slow
@pytest.mark.parametrize("scratch_buffer", [False, True], ids=["inmem", "scratch"])
def test_mf6model(scratch_buffer, function_tmpdir, targets, benchmark):
    test = TestFramework(
        name=simname,
        workspace=function_tmpdir,
        build=lambda t: build_models(t, scratch_buffer),
        check=lambda t: check_output(t),
        targets=targets,
        compare=None,
    )
    benchmark.pedantic(test.run, rounds=1, iterations=1)

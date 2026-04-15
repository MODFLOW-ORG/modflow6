"""
Test PRT-PRT exchange between two particle tracking models,
where one model is quad-refined. This exercises logic that
that transfers particles across the exchange boundary, and
verifies that particles are sent to the correct cells and
positions in the refined destination model.

Each base model domain has 2 cells, which together form a
4-cell line. The right model is refined with one level of
quad refinement, so it has 8 cells. Steady flow is set up
from left to right with CHD packages in the left-most and
right-most cells.

The flow and tracking models are connected on their shared
faces, with the right-most cell in the left model connected
to both left-most cells in the right model. IFLOWFACE is
configured in the CHD cells and in the exchange cells, on
the outer faces of the CHD cells and on the shared exchange
faces at the exchange boundary.

Particles are released from the left-most cell in the left
model. They should reach and cross the exchange boundary,
and terminate at the rightmost edge of the right model.
Two particles are released at Y coordinates such that one
will pass through the upper row of cells in the right model
and one through the lower row.

                      left  |  right
____________________________|____________________________
|             |             |      |      |      |      |
|      o------|-------------|------|------|------|----->x
|             |             |______|______|______|______|
|             |             |      |      |      |      |
|      o------|-------------|------|------|------|----->x
|_____________|_____________|______|______|______|______|
       o release
       x termination

"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridgen import Gridgen
from framework import TestFramework
from test_prt_prt_exg import get_model_name, plot_output

nlay = 1
nrow = 1
ncol = 2

delr = 1.0
delc = 1.0
delz = 1.0
top = 1.0
botm = [0.0]

head_left = 2.0
head_right = -1.0

hk = 1.0

porosity = 0.1
nper = 1
perlen = 10.0
nstp = 10

simname = "prtprtqref"
cases = [simname]


def get_refined_gridprops(test):
    """Create refined DISV grid for right model using gridgen."""
    workspace = test.workspace
    targets = test.targets

    # Create base grid for right model (1x1x2)
    ms = flopy.modflow.Modflow()
    dis = flopy.modflow.ModflowDis(
        ms,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # Create gridgen workspace
    gridgen_ws = workspace / "gridgen"
    gridgen_ws.mkdir(parents=True, exist_ok=True)

    # Create Gridgen object
    g = Gridgen(
        ms.modelgrid,
        model_ws=gridgen_ws,
        exe_name=targets["gridgen"],
    )

    # Refine entire right model with 1 level of refinement
    # This will turn 2 cells into 8 cells (each refined 2x2)
    polygon = [[(0.0, 0.0), (0.0, delc), (2 * delr, delc), (2 * delr, 0.0), (0.0, 0.0)]]
    g.add_refinement_features([polygon], "polygon", 1, range(nlay))
    g.build(verbose=False)

    return g.get_gridprops_disv()


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    mf6 = test.targets["mf6"]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        sim_ws=ws,
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(perlen, nstp, 1.0)],
    )
    flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        complexity="simple",
        outer_dvclose=1e-8,
        inner_dvclose=1e-8,
    )

    # Left GWF model - regular DIS grid (2 cells)
    gwfl = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfl"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwfl,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwfl, strt=head_left)
    flopy.mf6.ModflowGwfnpf(
        gwfl,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    flopy.mf6.ModflowGwfchd(
        gwfl,
        stress_period_data=[[(0, 0, 0), head_left, 1]],
        pname="chd_left",
        auxiliary=["IFLOWFACE"],
    )
    flopy.mf6.ModflowGwfoc(
        gwfl,
        head_filerecord=f"{get_model_name(idx, 'gwfl')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfl')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # Right GWF model - refined DISV grid (8 cells after refinement)
    gridprops = get_refined_gridprops(test)

    gwfr = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfr"),
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdisv(gwfr, **gridprops)
    flopy.mf6.ModflowGwfic(gwfr, strt=head_right)
    flopy.mf6.ModflowGwfnpf(
        gwfr,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    # CHD in rightmost cells of refined grid
    # From gridgen output: cells 6 and 8 (1-based) are at x=1.75-2.0 (rightmost)
    # In FloPy 0-based indexing: cells 5 and 7
    # IFLOWFACE 2 is east face for these cells
    flopy.mf6.ModflowGwfchd(
        gwfr,
        stress_period_data=[[(0, 5), head_right, 2], [(0, 7), head_right, 2]],
        pname="chd_right",
        auxiliary=["IFLOWFACE"],
    )
    flopy.mf6.ModflowGwfoc(
        gwfr,
        head_filerecord=f"{get_model_name(idx, 'gwfr')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfr')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    # GWF-GWF exchange
    # Left model cell (0,0,1) connects to right model cells 1 and 3 (0-based: 0 and 2)
    # These are the leftmost cells of the refined grid
    gwfgwf_data = [
        # Left cell (0,0,1) to refined right upper cell (cell 1, 0-based: 0)
        (
            (0, 0, ncol - 1),
            (0, 0),
            1,
            delr / 2.0,
            delr / 4.0,
            delc / 2.0,
            0.0,
            delr / 2.0,
        ),
        # Left cell (0,0,1) to refined right lower cell (cell 3, 0-based: 2)
        (
            (0, 0, ncol - 1),
            (0, 2),
            1,
            delr / 2.0,
            delr / 4.0,
            delc / 2.0,
            0.0,
            delr / 2.0,
        ),
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=get_model_name(idx, "gwfl"),
        exgmnameb=get_model_name(idx, "gwfr"),
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        filename="gwf-gwf.exg",
    )

    # Left PRT model
    prtl = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtl"))
    flopy.mf6.ModflowPrtdis(
        prtl,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowPrtmip(prtl, pname="mip", porosity=porosity)

    # Release 2 particles: one in upper half, one in lower half
    # This ensures particles go through different rows in refined right model
    releasepts = [
        (0, (0, 0, 0), 0.5, 0.75, 0.5),  # Upper particle
        (1, (0, 0, 0), 0.5, 0.25, 0.5),  # Lower particle
    ]
    flopy.mf6.ModflowPrtprp(
        prtl,
        pname="prp1",
        filename=f"{get_model_name(idx, 'prtl')}.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        extend_tracking=True,
        stoptime=0.4,
    )
    flopy.mf6.ModflowPrtoc(
        prtl,
        pname="oc",
        trackcsv_filerecord=f"{get_model_name(idx, 'prtl')}.trk.csv",
        dev_dump_event_trace=True,
    )
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(idx, "gwfl"),
        exgmnameb=get_model_name(idx, "prtl"),
        filename=f"{get_model_name(idx, 'gwfl')}.gwfprt",
    )

    # Right PRT model - refined DISV grid
    prtr = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtr"))
    flopy.mf6.ModflowPrtdisv(prtr, **gridprops)
    flopy.mf6.ModflowPrtmip(prtr, pname="mip", porosity=porosity)
    # Right model has no PRP, it just receives particles from left
    flopy.mf6.ModflowPrtoc(
        prtr,
        pname="oc",
        trackcsv_filerecord=f"{get_model_name(idx, 'prtr')}.trk.csv",
        dev_dump_event_trace=True,
    )
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(idx, "gwfr"),
        exgmnameb=get_model_name(idx, "prtr"),
        filename=f"{get_model_name(idx, 'gwfr')}.gwfprt",
    )

    # PRT-PRT exchange
    # For DISV grids, IFLOWFACE proceeds clockwise from first vertex
    # Left model DIS: face 3 is east (right) face
    # Right model DISV cells 1 and 3 both have west face = face 4
    # (verified from vertex coordinates - face 4 connects leftmost vertices)
    prtprt_data = [
        # Left cell (0,0,1) to refined upper right cell (cell 1, 0-based: 0)
        ((0, 0, ncol - 1), (0, 0), 1, 3, 4),  # DIS face 3 (east) to DISV face 4 (west)
        # Left cell (0,0,1) to refined lower right cell (cell 3, 0-based: 2)
        ((0, 0, ncol - 1), (0, 2), 1, 3, 4),  # DIS face 3 (east) to DISV face 4 (west)
    ]
    flopy.mf6.ModflowPrtprt(
        sim,
        exgtype="PRT6-PRT6",
        nexg=len(prtprt_data),
        exgmnamea=get_model_name(idx, "prtl"),
        exgmnameb=get_model_name(idx, "prtr"),
        gwfmodelname1=get_model_name(idx, "gwfl"),
        gwfmodelname2=get_model_name(idx, "gwfr"),
        exchangedata=prtprt_data,
        filename="prt-prt.exg",
        auxiliary=["IFLOWFACE1", "IFLOWFACE2"],
    )

    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename="prt.ems",
    )
    sim.register_solution_package(
        ems, [get_model_name(idx, "prtl"), get_model_name(idx, "prtr")]
    )

    return sim


def check_output(idx, test):
    ws = Path(test.workspace)
    hds_l = flopy.utils.HeadFile(ws / f"{get_model_name(idx, 'gwfl')}.hds")
    hds_r = flopy.utils.HeadFile(ws / f"{get_model_name(idx, 'gwfr')}.hds")
    head_l = hds_l.get_data().squeeze()
    head_r = hds_r.get_data().squeeze()
    assert np.all(head_l <= head_left)
    assert np.all(head_r >= head_right)

    pls_l = pd.read_csv(ws / f"{get_model_name(idx, 'prtl')}.trk.csv")
    pls_r = pd.read_csv(ws / f"{get_model_name(idx, 'prtr')}.trk.csv")
    assert len(pls_l) > 0
    assert len(pls_r) > 0

    irpts_l = pls_l["irpt"].unique()
    irpts_r = pls_r["irpt"].unique()
    assert set(irpts_l) == set(irpts_r)

    terminations = pls_r[pls_r.ireason == 3]
    assert len(terminations) == len(irpts_l)
    assert all(terminations.istatus == 2)

    final_x = terminations["x"]
    final_y = terminations["y"]
    final_z = terminations["z"]
    assert np.allclose(final_x, 2.0)
    assert all(np.isclose(y, 0.25) or np.isclose(y, 0.75) for y in final_y)
    assert np.allclose(final_z, 0.5)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        compare=None,
    )
    test.run()

"""
Test PRT-PRT exchange between two particle tracking models.

Each model domain has 2 cells, together forming a 4-cell line.
Steady flow is configured from left to right with CHD packages
in the left-most and right-most cells.

The flow and tracking models are connected on their shared face.
IFLOWFACE is configured in the CHD cells and in the exchange cells,
on the outer faces of the CHD cells and on the shared exchange face
at the exchange boundary.

A single particle is released from the left-most cell in the left
model. It should reach and cross the exchange and terminate at the
rightmost edge of the right model.

   left | right
________|________
| o-|---|---|-->x
|___|___|___|___|
  o release
  x termination

There are two test cases, one using a DIS grid, another with an
equivalent DISV grid.

"""

from pathlib import Path

import flopy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.gridutil import get_disv_kwargs
from framework import TestFramework

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

simname = "prtprt"
cases = [f"{simname}dis", f"{simname}disv"]


def get_model_name(idx, mdl):
    return f"{cases[idx]}_{mdl}"


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

    gwfl = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfl"),
        save_flows=True,
    )
    if "disv" in name:
        kwargs = get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
        flopy.mf6.ModflowGwfdisv(gwfl, **kwargs)
    else:
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
        stress_period_data=[[(0, 0), head_left, 4]]
        if "disv" in name
        else [[(0, 0, 0), head_left, 1]],
        pname="chd_left",
        auxiliary=["IFLOWFACE"],
    )
    flopy.mf6.ModflowGwfoc(
        gwfl,
        head_filerecord=f"{get_model_name(idx, 'gwfl')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfl')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwfr = flopy.mf6.ModflowGwf(
        sim,
        modelname=get_model_name(idx, "gwfr"),
        save_flows=True,
    )
    if "disv" in name:
        flopy.mf6.ModflowGwfdisv(
            gwfr, **get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
        )
    else:
        flopy.mf6.ModflowGwfdis(
            gwfr,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
    flopy.mf6.ModflowGwfic(gwfr, strt=head_right)
    flopy.mf6.ModflowGwfnpf(
        gwfr,
        save_specific_discharge=True,
        icelltype=0,
        k=hk,
    )
    flopy.mf6.ModflowGwfchd(
        gwfr,
        stress_period_data=[[(0, 1), head_right, 2]]
        if "disv" in name
        else [[(0, 0, 1), head_right, 3]],
        pname="chd_right",
        auxiliary=["IFLOWFACE"],
    )
    flopy.mf6.ModflowGwfoc(
        gwfr,
        head_filerecord=f"{get_model_name(idx, 'gwfr')}.hds",
        budget_filerecord=f"{get_model_name(idx, 'gwfr')}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    if "disv" in name:
        gwfgwf_data = [
            (
                (0, ncol - 1),
                (0, 0),
                1,
                delr / 2.0,
                delr / 2.0,
                delc,
                0.0,
                delr,
            )
        ]
    else:
        gwfgwf_data = [
            (
                (0, 0, ncol - 1),
                (0, 0, 0),
                1,
                delr / 2.0,
                delr / 2.0,
                delc,
                0.0,
                delr,
            )
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

    prtl = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtl"))
    if "disv" in name:
        flopy.mf6.ModflowPrtdisv(
            prtl, **get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
        )
        releasepts = [(0, (0, 0), 0.5, 0.5, 0.5)]
    else:
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
        releasepts = [(0, (0, 0, 0), 0.5, 0.5, 0.5)]
    flopy.mf6.ModflowPrtmip(prtl, pname="mip", porosity=porosity)
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

    prtr = flopy.mf6.ModflowPrt(sim, modelname=get_model_name(idx, "prtr"))
    if "disv" in name:
        flopy.mf6.ModflowPrtdisv(
            prtr, **get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, botm)
        )
    else:
        flopy.mf6.ModflowPrtdis(
            prtr,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
        )
    flopy.mf6.ModflowPrtmip(prtr, pname="mip", porosity=porosity)
    # right model has no PRP, it just receives particles from left
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

    if "disv" in name:
        prtprt_data = [
            # nodem1, nodem2, ihc, iflowface1, iflowface2
            ((0, ncol - 1), (0, 0), 1, 2, 4)
        ]
    else:
        prtprt_data = [
            # nodem1, nodem2, ihc, iflowface1, iflowface2
            ((0, 0, ncol - 1), (0, 0, 0), 1, 3, 1)
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

    final_x = terminations["x"].item()
    final_y = terminations["y"].item()
    final_z = terminations["z"].item()
    assert np.isclose(final_x, 2.0)
    assert np.isclose(final_y, 0.5)
    assert np.isclose(final_z, 0.5)


def plot_output(idx, test):
    name = test.name
    ws = test.workspace
    sim = test.sims[0]
    gwf_l_name = get_model_name(idx, "gwfl")
    gwf_r_name = get_model_name(idx, "gwfr")
    prt_l_name = get_model_name(idx, "prtl")
    prt_r_name = get_model_name(idx, "prtr")
    gwf_l = sim.get_model(gwf_l_name)
    gwf_r = sim.get_model(gwf_r_name)
    mg_l = gwf_l.modelgrid
    mg_r = gwf_r.modelgrid
    gwf_l_head_file = f"{gwf_l_name}.hds"
    gwf_r_head_file = f"{gwf_r_name}.hds"
    prt_l_trk_file = f"{prt_l_name}.trk.csv"
    prt_r_trk_file = f"{prt_r_name}.trk.csv"
    hds_l = flopy.utils.HeadFile(ws / gwf_l_head_file).get_data()
    hds_r = flopy.utils.HeadFile(ws / gwf_r_head_file).get_data()
    bud_l = gwf_l.output.budget()
    bud_r = gwf_r.output.budget()
    spdis_l = bud_l.get_data(text="DATA-SPDIS")[0]
    spdis_r = bud_r.get_data(text="DATA-SPDIS")[0]
    qs_l = flopy.utils.postprocessing.get_specific_discharge(spdis_l, gwf_l)
    qs_r = flopy.utils.postprocessing.get_specific_discharge(spdis_r, gwf_r)
    pls_l = pd.read_csv(ws / prt_l_trk_file, na_filter=False)
    pls_r = pd.read_csv(ws / prt_r_trk_file, na_filter=False)

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # left model
    pmv = flopy.plot.PlotMapView(modelgrid=mg_l, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds_l[0], alpha=0.1)
    pmv.plot_vector(qs_l[0], qs_l[1], normalize=True, color="white")
    plines_l = pls_l.groupby(["irpt"])
    for ipl, (irpt, pl) in enumerate(plines_l):
        pl.plot(
            title="left",
            kind="line",
            x="x",
            y="y",
            ax=ax[0],
            legend=False,
            color=cm.plasma(ipl / len(plines_l)),
        )

    # right model
    pmv = flopy.plot.PlotMapView(modelgrid=mg_r, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(hds_r[0], alpha=0.1)
    pmv.plot_vector(qs_r[0], qs_r[1], normalize=True, color="white")
    plines_r = pls_r.groupby(["irpt"])
    for ipl, (irpt, pl) in enumerate(plines_r):
        pl.plot(
            title="right",
            kind="line",
            x="x",
            y="y",
            ax=ax[1],
            legend=False,
            color=cm.plasma(ipl / len(plines_r)),
        )

    # view/save plot
    plt.show()
    plt.savefig(ws / f"{name}.png")


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

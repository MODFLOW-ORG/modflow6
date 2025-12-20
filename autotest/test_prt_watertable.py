"""
Tests particle tracking near the water table. This is a simplified
version of the Keating problem without the central aquitard. Those
particles dropped to the left of recharge zone should move left to
a stagnation point, along the surface of the water. Those released
within or to the right of the recharge zone should be moved by the
flow to terminate finally at the right-side constant head boundary.
"""

from os import environ
from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from modflow_devtools.misc import is_in_ci
from prt_test_utils import get_model_name

simname = "prtwt"
cases = [simname]

figure_size = (7.5, 3)
nlay = 80  # Number of layers
nrow = 1  # Number of rows
ncol = 400  # Number of columns
delr = 25.0  # Column width ($m$)
delc = 1.0  # Row width ($m$)
delz = 25.0  # Layer thickness ($m$)
top = 2000.0  # Top of model domain ($m$)
bottom = 0.0  # Bottom of model domain ($m$)
hka = 1.0e-12  # Permeability of aquifer ($m^2$)
h1 = 800.0  # Head on left side ($m$)
h2 = 100.0  # Head on right side ($m$)
recharge = 0.5  # Recharge ($kg/s$)
recharge_conc = 1.0  # Normalized recharge concentration (unitless)
alpha_l = 1.0  # Longitudinal dispersivity ($m$)
alpha_th = 1.0  # Transverse horizontal dispersivity ($m$)
alpha_tv = 1.0  # Transverse vertical dispersivity ($m$)
period1 = 730  # Length of first simulation period ($d$)
period2 = 29270.0  # Length of second simulation period ($d$)
porosity = 0.1  # Porosity of mobile domain (unitless)
seconds_to_days = 24.0 * 60.0 * 60.0
permeability_to_conductivity = 1000.0 * 9.81 / 1.0e-3 * seconds_to_days
hka = hka * permeability_to_conductivity
botm = [top - (k + 1) * delz for k in range(nlay)]
x = np.arange(0, 10000.0, delr) + delr / 2.0
hydraulic_conductivity = np.ones((nlay, nrow, ncol), dtype=float) * hka
rcol = []
for jcol in range(ncol):
    if 4200.0 <= x[jcol] <= 4800.0:
        rcol.append(jcol)
number_recharge_cells = len(rcol)
rrate = recharge * seconds_to_days / 1000.0
cell_area = delr * delc
rrate = rrate / (float(number_recharge_cells) * cell_area)
rchspd = {}
rchspd[0] = [[(0, 0, j), rrate, recharge_conc] for j in rcol]
rchspd[1] = [[(0, 0, j), rrate, 0.0] for j in rcol]


def build_gwf_sim(name, gwf_ws, mf6):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=gwf_ws, exe_name=mf6)
    tdis_ds = ((period1, 1, 1.0), (period2, 1, 1.0))
    flopy.mf6.ModflowTdis(sim, nper=len(tdis_ds), perioddata=tdis_ds)
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        complexity="complex",
        no_ptcrecord="all",
        outer_dvclose=1.0e-4,
        outer_maximum=2000,
        under_relaxation="dbd",
        linear_acceleration="BICGSTAB",
        under_relaxation_theta=0.7,
        under_relaxation_kappa=0.08,
        under_relaxation_gamma=0.05,
        under_relaxation_momentum=0.0,
        backtracking_number=20,
        backtracking_tolerance=2.0,
        backtracking_reduction_factor=0.2,
        backtracking_residual_limit=5.0e-4,
        inner_dvclose=1.0e-5,
        rcloserecord="0.0001 relative_rclose",
        inner_maximum=100,
        relaxation_factor=0.0,
        number_orthogonalizations=2,
        preconditioner_levels=8,
        preconditioner_drop_tolerance=0.001,
    )
    gwf_name = get_model_name(name, "gwf")
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwf_name, save_flows=True, newtonoptions=["newton"]
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_saturation=True,
        save_flows=True,
        icelltype=1,
        k=hydraulic_conductivity,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=600.0)
    chdspd = [[(k, 0, 0), h1] for k in range(nlay) if botm[k] < h1]
    chdspd += [[(k, 0, ncol - 1), h2] for k in range(nlay) if botm[k] < h2]
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspd,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfrch(
        gwf,
        stress_period_data=rchspd,
        auxiliary=["concentration"],
        pname="RCH-1",
    )

    head_filerecord = f"{gwf_name}.hds"
    budget_filerecord = f"{gwf_name}.bud"
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[
            ("HEAD", "ALL"),
            ("BUDGET", "ALL"),
        ],
    )
    return sim


def build_prt_sim(name, gwf, prt_ws, mf6):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=prt_ws,
        exe_name=mf6,
        continue_=True,
    )
    tdis_ds = ((period1, 1, 1.0), (period2, 1, 1.0))
    flopy.mf6.ModflowTdis(sim, nper=len(tdis_ds), perioddata=tdis_ds)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    dis = flopy.mf6.ModflowGwtdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)
    nns = range(prt.modelgrid.ncpl)
    ids = prt.modelgrid.get_lrc(nns)
    ccs = list(
        zip(prt.modelgrid.xcellcenters.ravel(), prt.modelgrid.ycellcenters.ravel())
    )
    prpdata = [(nn, *ids[nn], *ccs[nn], 0.5) for nn in nns]
    prp = flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1a",
        filename=f"{prt_name}_1a.prp",
        nreleasepts=len(prpdata),
        packagedata=prpdata,
        releasetimes=[(0,)],
        nreleasetimes=1,
        exit_solve_tolerance=1e-5,
        extend_tracking=False,
        local_z=True,
    )
    budget_record = [f"{prt_name}.cbc"]
    track_record = [f"{prt_name}.trk"]
    trackcsv_record = [f"{prt_name}.trk.csv"]
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        budget_filerecord=budget_record,
        track_filerecord=track_record,
        trackcsv_filerecord=trackcsv_record,
        saverecord=[("BUDGET", "ALL")],
    )
    pd = [
        ("GWFGRID", Path(f"../{Path(gwf.model_ws).name}/{gwf_name}.dis.grb")),
        ("GWFHEAD", Path(f"../{Path(gwf.model_ws).name}/{gwf_name}.hds")),
        ("GWFBUDGET", Path(f"../{Path(gwf.model_ws).name}/{gwf_name}.bud")),
    ]
    fmi = flopy.mf6.ModflowPrtfmi(prt, packagedata=pd)
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt_name])
    return sim


def build_models(idx, test):
    """Build both GWF and PRT simulations for the given test case."""

    gwf_sim = build_gwf_sim(test.name, test.workspace / "gwf", test.targets["mf6"])
    prt_sim = build_prt_sim(
        test.name,
        gwf_sim.get_model(),
        test.workspace / "prt",
        test.targets["mf6"],
    )
    return gwf_sim, prt_sim


def check_output(idx, test, snapshot):
    name = test.name
    prt_ws = test.workspace / "prt"
    prt_name = get_model_name(name, "prt")

    trackcsv_file = f"{prt_name}.trk.csv"
    trackcsv_path = prt_ws / trackcsv_file
    pls = pd.read_csv(trackcsv_path)

    strtpts = pls[pls.ireason == 0]
    endpts = pls[pls.ireason == 3]
    n_particles = len(strtpts)
    assert len(endpts) == n_particles

    if snapshot is not None:
        actual = (
            pls.drop(["name", "icell"], axis=1, errors="ignore")
            .round(1)
            .reset_index(drop=True)
        )
        assert snapshot == actual.to_records(index=False)


def plot_output(idx, test):
    sim_gwf, sim_prt = test.sims
    gwf_name = get_model_name(test.name, "gwf")
    prt_name = get_model_name(test.name, "prt")
    gwf = sim_gwf.get_model(gwf_name)
    fig, ax = plt.subplots(1, 1, figsize=figure_size, dpi=300, tight_layout=True)
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    pxs.plot_ibound()
    pxs.plot_bc(ftype="RCH", color="red")
    pxs.plot_bc(ftype="CHD")
    ax.set_aspect(1.0)
    pathlines = pd.read_csv(sim_prt.sim_path / f"{prt_name}.trk.csv")
    for _, pl in pathlines.groupby("irpt"):
        pl.plot("x", "z", lw=0.2, ax=ax, alpha=0.4, legend=False)
    plt.show()


@pytest.mark.snapshot
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, array_snapshot, plot):
    # skip snapshot in CI unless gfortran. flow model results vary for
    # intel vs gcc compilers when built with optimizations
    if is_in_ci() and "gfortran" not in environ.get("FC", "").lower():
        array_snapshot = None
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t, array_snapshot),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()

"""
Tests particle mass budget tracking.
"""

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from flopy.utils.binaryfile import HeadFile
from framework import TestFramework
from prt_test_utils import (
    get_model_name,
)

simname = "prtbud"
cases = [simname]

# model names
gwf_name = get_model_name(simname, "gwf")
gwt_name = get_model_name(simname, "gwt")
prt_name = get_model_name(simname, "prt")

# tdis data
years = 550.0
year_delt = 365.25
perlen = year_delt * years
nstp = 10
tsmult = 1.0
tdis_spd = [[perlen, nstp, tsmult]]
nper = len(tdis_spd)

# grid data
nlay, nrow, ncol = 5, 101, 101
delr, delc = 100.0, 100.0
top = 1.0
botm = np.linspace(-10, -100, nlay)

# gwf chd data
chd_spd = [[(0, i, 0), 1.0, 6, -1] for i in range(nrow)]
chd_spd += [[(0, i, ncol - 1), 0.0, 6, -1] for i in range(nrow)]

# gwf maw data
maw_spd = [((2, 50, 50), -5000.0, 0, 0), ((4, 50, 50), -5000.0, 0, 0)]

# gwf output file names
gwf_budget_file = f"{gwf_name}.bud"
gwf_head_file = f"{gwf_name}.hds"

# gwt src data
gwt_srcs = []
for k in [0, 1, 2, 3]:  # range(nlay):
    gwt_srcs += [(k, i, 1) for i in range(nrow)]
src_rate = 4.9779105e-4
src_spd = []
for cid in gwt_srcs:
    src_spd.append([cid, "SRCRATE"])
src_spd = {0: src_spd}
src_tsdata = [
    (0.0, src_rate),
    (perlen / 100.0, 0.0),
    (perlen, 0.0),
]

# prt prp data
prt_nodes = []
for k in [0, 1, 2, 3]:
    prt_nodes += [(k, i, 1) for i in range(nrow)]
particle_data = flopy.modpath.ParticleData(
    prt_nodes,
    drape=0,
    structured=True,
)

# prt mip data
izone = np.zeros((nlay, nrow, ncol), dtype=int)
well_zone = (1, 2)
for idx, k in enumerate((2, 4)):
    izone[k, 50, 50] = well_zone[idx]

# prt oc data
tracktimes = list(range(0, 72000, 1000))

# prt output file names
prt_listfile = f"{prt_name}.lst"
prt_budgetfile = f"{prt_name}.cbb"
prt_trackfile = f"{prt_name}.trk"
prt_trackcsvfile = f"{prt_name}.trk.csv"


def build_gwf_sim(gwf_ws, mf6):
    sim = flopy.mf6.MFSimulation(sim_name=gwf_name, sim_ws=gwf_ws, exe_name=mf6)
    tdis = flopy.mf6.ModflowTdis(sim, nper=1, perioddata=tdis_spd)
    ims = flopy.mf6.ModflowIms(sim, linear_acceleration="bicgstab", complexity="simple")
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        save_flows=True,
        newtonoptions="newton under_relaxation",
    )
    dis = flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=nrow, ncol=ncol, delr=100, delc=100, top=top, botm=botm
    )
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=1, save_specific_discharge=True, save_saturation=True, k=100.0
    )
    sto = flopy.mf6.ModflowGwfsto(
        gwf, iconvert=1, steady_state={0: False}, transient={0: True}
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary=["iface", "iflowface"],
        maxbound=len(chd_spd),
        stress_period_data=chd_spd,
        pname="chd",
    )
    maw = flopy.mf6.ModflowGwfwel(
        gwf,
        auxiliary=["iface", "iflowface"],
        maxbound=len(maw_spd),
        stress_period_data=maw_spd,
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=gwf_head_file,
        budget_filerecord=gwf_budget_file,
        printrecord=[("budget", "all")],
        saverecord=[("head", "all"), ("budget", "all")],
    )
    return sim


def build_gwt_sim(gwf_ws, gwt_ws, mf6):
    sim = flopy.mf6.MFSimulation(sim_name=gwt_name, sim_ws=gwt_ws, exe_name=mf6)
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, nper=nper, perioddata=tdis_spd)
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwt_name, print_input=True)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0)
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt,
        xt3d_off=True,
        alh=250.0,
        ath1=25.0,
        ath2=25.0,
    )
    ssm = flopy.mf6.ModflowGwtssm(gwt)
    src = flopy.mf6.ModflowGwtsrc(
        gwt,
        stress_period_data=src_spd,
        timeseries={
            "timeseries": src_tsdata,
            "time_series_namerecord": "SRCRATE",
            "interpolation_methodrecord": "STEPWISE",
        },
    )
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwt_name}.cbb",
        concentration_filerecord=f"{gwt_name}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "ALL")],
    )
    fmi = flopy.mf6.ModflowGwtfmi(
        gwt,
        packagedata=[
            ("GWFHEAD", gwf_ws / gwf_head_file),
            ("GWFBUDGET", gwf_ws / gwf_budget_file),
        ],
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        pname="ims",
        filename=f"{gwt_name}.ims",
        print_option="summary",
        inner_maximum=300,
        linear_acceleration="bicgstab",
        inner_dvclose=1e-9,
    )
    sim.register_solution_package(ims, [gwt.name])
    return sim


def build_prt_sim(gwf_ws, prt_ws, mf6):
    sim = flopy.mf6.MFSimulation(sim_name=prt_name, sim_ws=prt_ws, exe_name=mf6)
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, nper=nper, perioddata=tdis_spd)
    prt = flopy.mf6.ModflowPrt(
        sim, modelname=prt_name, print_input=True, save_flows=True
    )
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    mip = flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=0.1, izone=izone)
    releasepts = list(particle_data.to_prp(prt.modelgrid))
    prp = flopy.mf6.ModflowPrtprp(
        prt,
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: ["FIRST"]},
        exit_solve_tolerance=1e-5,
        extend_tracking=False,
        print_input=True,
    )
    oc = flopy.mf6.ModflowPrtoc(
        prt,
        budget_filerecord=[prt_budgetfile],
        track_filerecord=[prt_trackfile],
        trackcsv_filerecord=[prt_trackcsvfile],
        ntracktimes=len(tracktimes),
        tracktimes=[(t,) for t in tracktimes],
        printrecord=[("BUDGET", "ALL")],
        saverecord=[("BUDGET", "ALL")],
    )
    fmi = flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", gwf_ws / gwf_head_file),
            ("GWFBUDGET", gwf_ws / gwf_budget_file),
        ],
    )
    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(idx, test):
    gwf_sim = build_gwf_sim(test.workspace / "gwf", test.targets["mf6"])
    gwt_sim = build_gwt_sim(
        test.workspace / "gwf", test.workspace / "gwt", test.targets["mf6"]
    )
    prt_sim = build_prt_sim(
        test.workspace / "gwf", test.workspace / "prt", test.targets["mf6"]
    )
    return gwf_sim, gwt_sim, prt_sim


def check_output(idx, test):
    gwf_sim = test.sims[0]
    gwt_sim = test.sims[1]
    prt_sim = test.sims[2]
    gwf_ws = test.workspace / "gwf"
    gwt_ws = test.workspace / "gwt"
    prt_ws = test.workspace / "prt"
    gwf = gwf_sim.get_model(gwf_name)
    gwt = gwt_sim.get_model(gwt_name)
    prt = prt_sim.get_model(prt_name)
    grid = gwf.modelgrid

    # check cumulative budget from list file
    prt_lst = flopy.utils.Mf6ListBudget(
        prt_ws / prt_listfile, budgetkey="MASS BUDGET FOR ENTIRE MODEL"
    )
    for key in [
        "STORAGE_IN",
        "PRP_IN",
        "TERMINATION_IN",
        "STORAGE_OUT",
        "PRP_OUT",
        "TERMINATION_OUT",
    ]:
        assert key in prt_lst.get_record_names()
    prt_bud_cum = prt_lst.get_data()
    prp_in = next(
        iter([term[1] for term in prt_bud_cum if term[2].decode() == "PRP_IN"])
    )
    prp_out = next(
        iter([term[1] for term in prt_bud_cum if term[2].decode() == "PRP_OUT"])
    )
    term_in = next(
        iter([term[1] for term in prt_bud_cum if term[2].decode() == "TERMINATION_IN"])
    )
    term_out = next(
        iter([term[1] for term in prt_bud_cum if term[2].decode() == "TERMINATION_OUT"])
    )
    pct_dscr = next(
        iter(
            [
                term[1]
                for term in prt_bud_cum
                if term[2].decode() == "PERCENT_DISCREPANCY"
            ]
        )
    )
    assert np.isclose(prp_out, 0.0)
    assert np.isclose(term_in, 0.0)
    assert np.isclose(pct_dscr, 0.0)
    assert np.isclose(prp_in + term_out, 0.0)

    # check cell by cell budget from binary budget file
    prt_pls = pd.read_csv(prt_ws / prt_trackcsvfile, na_filter=False)
    prt_bud = flopy.utils.CellBudgetFile(prt_ws / prt_budgetfile, precision="double")
    prp_bud = prt_bud.get_data(text="PRP")
    assert len(prp_bud) == nstp
    assert len(prp_bud[0]["q"]) == prt_pls.irpt.unique().size
    # expect positive mass flow on release in first time step
    assert np.all(prp_bud[0]["q"] > 0)
    # and none in all subsequent steps
    for i in range(1, nstp):
        assert np.all(prp_bud[i]["q"] == 0)

    # TODO check storage and termination terms once added to binary budget file


def plot_output(idx, test):
    gwf_sim = test.sims[0]
    gwt_sim = test.sims[1]
    prt_sim = test.sims[2]
    gwf_ws = test.workspace / "gwf"
    gwt_ws = test.workspace / "gwt"
    prt_ws = test.workspace / "prt"
    gwf = gwf_sim.get_model(gwf_name)
    gwt = gwt_sim.get_model(gwt_name)
    prt = prt_sim.get_model(prt_name)
    grid = gwf.modelgrid

    # get head, budget, and spdis from GWF model
    hds = HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # get concentration from GWT model
    cobj = gwt.output.concentration()
    conc = cobj.get_alldata()
    pl = flopy.plot.PlotMapView(model=gwf)
    pl.ax.set_title("Concentration")
    pl.plot_grid(lw=0.5, alpha=0.5)
    pl.plot_array(hds[0], alpha=0.1)
    pl.contour_array(
        conc[99, 0],
        colors="blue",
        linestyles="-",
    )
    plt.show()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()

"""
Tests particle mass budget tracking.
"""

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from matplotlib import cm
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


def check_cumulative_prt_budget(path, nparticles):
    prt_lst = flopy.utils.Mf6ListBudget(path, budgetkey="MASS BUDGET FOR ENTIRE MODEL")

    expected_terms = [
        "PRP_IN",
        "PRP_OUT",
        "STORAGE_IN",
        "STORAGE_OUT",
        "TERMINATION_IN",
        "TERMINATION_OUT",
    ]
    actual_terms = prt_lst.get_record_names()
    for term in expected_terms:
        assert term in actual_terms

    prt_bud_cum = prt_lst.get_data()

    def get_budget_term(term_name):
        matches = [term[1] for term in prt_bud_cum if term[2].decode() == term_name]
        assert len(matches) == 1
        return matches[0]

    prp_in = get_budget_term("PRP_IN")
    prp_out = get_budget_term("PRP_OUT")
    sto_in = get_budget_term("STORAGE_IN")
    sto_out = get_budget_term("STORAGE_OUT")
    term_in = get_budget_term("TERMINATION_IN")
    term_out = get_budget_term("TERMINATION_OUT")
    pct_dscr = get_budget_term("PERCENT_DISCREPANCY")

    assert np.isclose(prp_in, nparticles)  # all particles released
    assert np.isclose(prp_out, 0.0)  # no mass out of prp term
    assert sto_in >= 0.0
    assert np.isclose(sto_in, -sto_out)  # storage budget balance
    assert np.isclose(term_in, 0.0)  # no mass into termination term
    assert np.isclose(term_out, -nparticles)  # all particles terminated
    assert np.isclose(pct_dscr, 0.0, atol=1e-6)  # overall budget balance
    assert np.isclose(
        prp_in + sto_in + term_in, -(prp_out + sto_out + term_out), rtol=1e-6
    )


def check_cell_by_cell_budget(path, nparticles):
    prt_bud = flopy.utils.CellBudgetFile(path, precision="double")
    prp_bud = prt_bud.get_data(text="PRP")
    sto_bud = prt_bud.get_data(text="STORAGE")
    trm_bud = prt_bud.get_data(text="TERMINATION")

    assert len(prp_bud) == nstp
    assert len(sto_bud) == nstp - 1
    assert len(trm_bud) == nstp - 1

    rls_bud = prp_bud[0]

    assert len(rls_bud) == nparticles  # correct number of particles released
    assert set(rls_bud["node2"]) == set(range(1, nparticles + 1))  # all released
    assert np.all(rls_bud["q"] > 0)  # all cell release rates positive (sparse data)
    # cell release rates in first time step
    for i, record in enumerate(prp_bud[0]):
        ic, ip, rate = record["node"], record["node2"], record["q"]
        assert ic > 0  # valid cell number
        assert ip == i + 1  # valid particle id
        assert np.isclose(rate, 1 / (perlen / nstp))
    # expect no further release
    for step in range(1, nstp):
        assert len(prp_bud[step]) == 0 or np.all(prp_bud[step]["q"] == 0)

    # 2nd time step, all particles in storage
    assert sum(sto_bud[0]["q"]) == nparticles

    # expect all particles terminated by end of simulation
    assert sum([sum(trm_bud[i]["q"]) for i in range(len(trm_bud))]) == nparticles


def check_output(idx, test):
    prt_ws = test.workspace / "prt"
    prt_pls = pd.read_csv(prt_ws / prt_trackcsvfile, na_filter=False)

    nparticles = prt_pls.irpt.unique().size
    assert len(prt_nodes) == nparticles

    check_cumulative_prt_budget(prt_ws / prt_listfile, nparticles)
    check_cell_by_cell_budget(prt_ws / prt_budgetfile, nparticles)


def plot_output(idx, test):
    name = test.name
    gwf_ws = test.workspace
    prt_ws = test.workspace / "prt"
    mp7_ws = test.workspace / "mp7"
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    mp7_name = get_model_name(name, "mp7")
    gwf_sim = test.sims[0]
    gwf = gwf_sim.get_model(gwf_name)
    mg = gwf.modelgrid

    # check mf6 output files exist
    gwf_head_file = f"{gwf_name}.hds"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    mp7_pathline_file = f"{mp7_name}.mppth"

    # load mf6 pathline results
    mf6_pls = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)

    # load mp7 pathline results
    plf = flopy.utils.PathlineFile(mp7_ws / mp7_pathline_file)
    mp7_pls = pd.DataFrame(
        plf.get_destination_pathline_data(range(mg.nnodes), to_recarray=True)
    )
    # convert zero-based to one-based indexing in mp7 results
    mp7_pls["particlegroup"] = mp7_pls["particlegroup"] + 1
    mp7_pls["node"] = mp7_pls["node"] + 1
    mp7_pls["k"] = mp7_pls["k"] + 1

    # extract head, budget, and specific discharge results from GWF model
    hds = flopy.utils.HeadFile(gwf_ws / gwf_head_file).get_data()
    bud = gwf.output.budget()
    spdis = bud.get_data(text="DATA-SPDIS")[0]
    qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

    # set up plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 10))
    for a in ax:
        a.set_aspect("equal")

    # plot mf6 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[0])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mf6_plines = mf6_pls.groupby(["iprp", "irpt", "trelease"])
    for ipl, ((iprp, irpt, trelease), pl) in enumerate(mf6_plines):
        pl.plot(
            title="MF6 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[0],
            legend=False,
            color=cm.plasma(ipl / len(mf6_plines)),
        )

    # plot mp7 pathlines in map view
    pmv = flopy.plot.PlotMapView(modelgrid=mg, ax=ax[1])
    pmv.plot_grid()
    pmv.plot_array(hds[0], alpha=0.1)
    pmv.plot_vector(qx, qy, normalize=True, color="white")
    mp7_plines = mp7_pls.groupby(["particleid"])
    for ipl, (pid, pl) in enumerate(mp7_plines):
        pl.plot(
            title="MP7 pathlines",
            kind="line",
            x="x",
            y="y",
            ax=ax[1],
            legend=False,
            color=cm.plasma(ipl / len(mp7_plines)),
        )

    # view/save plot
    plt.show()
    plt.savefig(prt_ws / f"test_{simname}.png")


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

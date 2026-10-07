"""
Tests the STRANDED_MASS option of the LKT and SFT Packages, which holds the
solute of a lake or reach that goes dry and returns it as the feature rewets.

Cases:
  - lake          : a lake on an impermeable bed evaporates dry and is refilled
                    by rainfall; all of its solute is held while it is dry and
                    is returned in proportion to the volume regained, relative
                    to its largest volume before it went dry.
  - lake_nooption : the same lake without the option cannot be solved once the
                    lake has no water (xfail).
  - reach         : the inflow to three reaches stops, so they go dry, and then
                    resumes with clean water; the solute of each reach is held
                    and returned when the reach refills.
  - energy        : the option is an error for the LKE Package (xfail).
"""

import re

import flopy
import numpy as np
from framework import TestFramework

cinit = 100.0


def add_ims(sim, fname):
    return flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="bicgstab",
        filename=fname,
    )


def add_gwf(sim, nlay, botm, idomain):
    gwf = flopy.mf6.ModflowGwf(sim, modelname="gwf", save_flows=True)
    sim.register_ims_package(add_ims(sim, "gwf.ims"), ["gwf"])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=1,
        ncol=3,
        delr=100.0,
        delc=100.0,
        top=10.0,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=4.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-5, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data=[[(nlay - 1, 0, 0), 4.0], [(nlay - 1, 0, 2), 4.0]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord="gwf.hds",
        budget_filerecord="gwf.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def add_gwt(sim, nlay, botm, idomain):
    gwt = flopy.mf6.ModflowGwt(sim, modelname="gwt", save_flows=True)
    sim.register_ims_package(add_ims(sim, "gwt.ims"), ["gwt"])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=1,
        ncol=3,
        delr=100.0,
        delc=100.0,
        top=10.0,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtadv(gwt)
    flopy.mf6.ModflowGwtssm(gwt)
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord="gwt.ucn",
        budget_filerecord="gwt.cbc",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )
    flopy.mf6.ModflowGwfgwt(sim, exgtype="GWF6-GWT6", exgmnamea="gwf", exgmnameb="gwt")
    return gwt


def lake_model(ws, strand):
    """A lake that evaporates dry over the second period and refills in the
    fourth, on an impermeable bed so that it loses water only to evaporation."""
    sim = flopy.mf6.MFSimulation(sim_name="lake", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=4, perioddata=[(10.0, 10, 1.0)] * 4)
    idomain = np.ones((2, 1, 3), dtype=int)
    idomain[0, 0, 1] = 0
    botm = [5.0, 0.0]
    gwf = add_gwf(sim, 2, botm, idomain)
    flopy.mf6.ModflowGwflak(
        gwf,
        pname="LAK-1",
        save_flows=True,
        nlakes=1,
        noutlets=0,
        packagedata=[[0, 6.0, 1]],
        connectiondata=[[0, 0, (1, 0, 1), "vertical", 0.0, 0.0, 0.0, 0.0, 0.0]],
        perioddata={
            0: [[0, "evaporation", 0.0]],
            1: [[0, "evaporation", 0.2]],
            2: [[0, "evaporation", 0.0]],
            3: [[0, "rainfall", 0.15]],
        },
        budget_filerecord="gwf.lak.bud",
    )
    gwt = add_gwt(sim, 2, botm, idomain)
    kwargs = {}
    if strand:
        kwargs.update(
            stranded_mass=True,
            observations={"lkt.obs.csv": [("held", "stranded", 1)]},
        )
    flopy.mf6.ModflowGwtlkt(
        gwt,
        flow_package_name="LAK-1",
        packagedata=[[0, cinit]],
        print_flows=True,
        budget_filerecord="gwt.lkt.bud",
        concentration_filerecord="gwt.lkt.conc",
        pname="LKT-1",
        **kwargs,
    )
    return sim


def reach_model(ws, strand):
    """Three reaches whose inflow stops in the second period, so that they go
    dry, and resumes with clean water in the third."""
    sim = flopy.mf6.MFSimulation(sim_name="reach", sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=3, perioddata=[(5.0, 5, 1.0)] * 3)
    idomain = np.ones((1, 1, 3), dtype=int)
    botm = [0.0]
    gwf = add_gwf(sim, 1, botm, idomain)
    # streambed conductance is zero, so the reaches exchange no water with the
    # aquifer
    pd = [
        [i, (0, 0, i), 100.0, 10.0, 1e-3, 9.0, 1.0, 0.0, 0.03, 1 if i != 1 else 2]
        + [1.0, 0]
        for i in range(3)
    ]
    flopy.mf6.ModflowGwfsfr(
        gwf,
        pname="SFR-1",
        save_flows=True,
        nreaches=3,
        packagedata=pd,
        connectiondata=[[0, -1], [1, 0, -2], [2, 1]],
        perioddata={
            0: [[0, "inflow", 100.0]],
            1: [[0, "inflow", 0.0]],
            2: [[0, "inflow", 100.0]],
        },
        budget_filerecord="gwf.sfr.bud",
    )
    gwt = add_gwt(sim, 1, botm, idomain)
    flopy.mf6.ModflowGwtsft(
        gwt,
        flow_package_name="SFR-1",
        packagedata=[[i, 0.0] for i in range(3)],
        reachperioddata={0: [[0, "inflow", cinit]], 2: [[0, "inflow", 0.0]]},
        print_flows=True,
        budget_filerecord="gwt.sft.bud",
        concentration_filerecord="gwt.sft.conc",
        pname="SFT-1",
        stranded_mass=strand,
    )
    return sim


def feature_series(ws, flowbud, transportbud, conc, strand):
    """Volume, concentration, and held mass of each feature at each time."""
    fb = flopy.utils.CellBudgetFile(ws / flowbud, precision="double")
    vol = np.array([r["VOLUME"] for r in fb.get_data(text="STORAGE")])
    cf = flopy.utils.HeadFile(ws / conc, text="CONCENTRATION")
    c = np.array([cf.get_data(totim=t).flatten() for t in cf.get_times()])
    held = np.zeros_like(vol)
    if strand:
        tb = flopy.utils.CellBudgetFile(ws / transportbud, precision="double")
        held = np.array([r["MASS"] for r in tb.get_data(text="STRANDED")])
    return vol, c, held


def assert_package_budget_closes(ws, pname):
    text = (ws / "gwt.lst").read_text()
    disc = re.findall(
        rf"{pname} BUDGET.*?PERCENT DISCREPANCY\s*=\s*([-+0-9.Ee]+)", text, re.S
    )
    assert disc, f"no {pname} budget in gwt.lst"
    disc = np.abs(np.array(disc, dtype=float))
    assert np.all(disc < 1e-2), f"{pname} budget does not close {disc}"


def test_lake(function_tmpdir, targets):
    def build(test):
        return lake_model(test.workspace, True)

    def check(test):
        ws = test.workspace
        assert_package_budget_closes(ws, "LKT-1")
        vol, c, held = feature_series(
            ws, "gwf.lak.bud", "gwt.lkt.bud", "gwt.lkt.conc", True
        )
        vol, c, held = vol[:, 0], c[:, 0], held[:, 0]
        total = vol * c + held
        mass0 = 10000.0 * cinit
        # the lake and the mass it holds keep all of the solute
        assert np.allclose(total, mass0, rtol=1e-6), f"mass not conserved {total}"
        # all of it is held while the lake is dry
        dry = vol <= 0.0
        assert dry.any(), "the lake did not go dry"
        assert np.allclose(held[dry], mass0, rtol=1e-6), f"held {held[dry]}"
        # the first refill step regains 1,500 of the 10,000 m3 the lake held
        # before it went dry, its largest volume, so 15 percent of it returns
        first = np.argmax(~dry & (np.arange(dry.size) > np.argmax(dry)))
        assert np.isclose(vol[first], 1500.0), f"volume {vol[first]}"
        assert np.isclose(held[first], 0.85 * mass0, rtol=1e-6), f"{held[first]}"
        # the rest has returned once the lake regains 10,000 m3
        full = first + np.argmax(vol[first:] >= 10000.0)
        assert np.isclose(held[full - 1], 0.1 * mass0, rtol=1e-6), (
            f"held before the lake refills {held[full - 1]}"
        )
        assert np.isclose(held[full], 0.0, atol=1e-6), f"held {held[full]}"
        # the stranded observation is the rate to the held mass, so over the
        # 1-d time step in which the lake dries it is the whole mass
        obs = np.genfromtxt(ws / "lkt.obs.csv", delimiter=",", names=True)
        assert np.isclose(obs["HELD"].min(), -mass0 / 1.0, rtol=1e-6), (
            f"stranded observation {obs['HELD'].min()} expected {-mass0}"
        )

    TestFramework(
        name="lake",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


def test_lake_nooption(function_tmpdir, targets):
    def build(test):
        return lake_model(test.workspace, False)

    def check(test):
        assert "convergence failure" in (test.workspace / "mfsim.lst").read_text()

    TestFramework(
        name="lake",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        xfail=True,
    ).run()


def test_reach(function_tmpdir, targets):
    def build(test):
        return reach_model(test.workspace, True)

    def check(test):
        ws = test.workspace
        assert_package_budget_closes(ws, "SFT-1")
        vol, c, held = feature_series(
            ws, "gwf.sfr.bud", "gwt.sft.bud", "gwt.sft.conc", True
        )
        # the reaches go dry in the first time step of the second period and
        # hold the solute they had at the end of the first
        last_wet, first_dry = 4, 5
        assert np.all(vol[first_dry] == 0.0), f"reaches not dry {vol[first_dry]}"
        expected = vol[last_wet] * c[last_wet]
        assert np.allclose(held[first_dry], expected, rtol=1e-6), (
            f"held {held[first_dry]} expected {expected}"
        )
        # each reach refills to its largest volume before it went dry, so all
        # of it returns in the first time step of the third period
        assert np.allclose(held[10], 0.0, atol=1e-6), f"held after refill {held[10]}"

    TestFramework(
        name="reach",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
    ).run()


def test_energy(function_tmpdir, targets):
    def build(test):
        ws = test.workspace
        sim = flopy.mf6.MFSimulation(sim_name="energy", sim_ws=ws, exe_name="mf6")
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
        idomain = np.ones((2, 1, 3), dtype=int)
        idomain[0, 0, 1] = 0
        botm = [5.0, 0.0]
        gwf = add_gwf(sim, 2, botm, idomain)
        flopy.mf6.ModflowGwflak(
            gwf,
            pname="LAK-1",
            nlakes=1,
            noutlets=0,
            packagedata=[[0, 6.0, 1]],
            connectiondata=[[0, 0, (1, 0, 1), "vertical", 0.0, 0.0, 0.0, 0.0, 0.0]],
        )
        gwe = flopy.mf6.ModflowGwe(sim, modelname="gwe")
        sim.register_ims_package(add_ims(sim, "gwe.ims"), ["gwe"])
        flopy.mf6.ModflowGwedis(
            gwe,
            nlay=2,
            nrow=1,
            ncol=3,
            delr=100.0,
            delc=100.0,
            top=10.0,
            botm=botm,
            idomain=idomain,
        )
        flopy.mf6.ModflowGweic(gwe, strt=10.0)
        flopy.mf6.ModflowGweest(
            gwe, porosity=0.3, heat_capacity_solid=800.0, density_solid=2600.0
        )
        flopy.mf6.ModflowGweadv(gwe)
        flopy.mf6.ModflowGwessm(gwe)
        flopy.mf6.ModflowGwelke(
            gwe,
            flow_package_name="LAK-1",
            packagedata=[[0, 10.0, 1.0, 1.0]],
            pname="LKE-1",
            filename="gwe.lke",
        )
        flopy.mf6.ModflowGwfgwe(
            sim, exgtype="GWF6-GWE6", exgmnamea="gwf", exgmnameb="gwe"
        )
        sim.write_simulation(silent=True)
        # the option is not part of the LKE input definition, so it is added
        # to the file directly
        fpth = ws / "gwe.lke"
        text = fpth.read_text()
        fpth.write_text(
            text.replace("BEGIN options", "BEGIN options\n  STRANDED_MASS", 1)
        )
        return sim

    def check(test):
        listing = " ".join((test.workspace / "mfsim.lst").read_text().split())
        assert "STRANDED_MASS is not supported for energy transport" in listing

    TestFramework(
        name="energy",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check,
        xfail=True,
        overwrite=False,
    ).run()

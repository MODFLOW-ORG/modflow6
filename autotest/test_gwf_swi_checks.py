"""
Regression tests for SWI input checks and for SWI combined with other GWF
features.

swi-tvs
  transient single-fluid strip whose STO has a TVS file that changes SY in
  the second stress period; the run must finish and conserve mass (SWI used
  to dereference the previous-SS and previous-SY arrays without setting
  them).
swi-tr2ss
  transient period followed by a steady-state period; the SWI storage term
  must be zero in the steady-state period rather than repeat the last
  transient value.
swi-tvg-start
  saltwater head supplied through TVG at a value consistent with the
  starting freshwater head; the interface does not move, so the SWI storage
  term must be zero from the first time step (the first TVG value used to be
  ignored as the initial saltwater head).
swi-picard
  SWI without Newton in a fully fresh confined strip; the interface is below
  the cells, so the flow must equal the flow without SWI.
swi-xt3d, swi-hfb
  SWI with a package it does not support; the run must stop with an error.
swi2f-grid, swi2f-sy
  two-fluid models whose grids (cell bottoms) or specific yields differ; the
  run must stop with an error.
swi2f-tvg
  two-fluid freshwater model with a TVG file for the saltwater head, which is
  simulated in two-fluid mode; the run must stop with an error.
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework

cases = [
    "swi-tvs",
    "swi-tr2ss",
    "swi-tvg-start",
    "swi-picard",
    "swi-xt3d",
    "swi-hfb",
    "swi2f-grid",
    "swi2f-sy",
    "swi2f-tvg",
]
xfail = [False, False, False, False, True, True, True, True, True]
errors = {
    "swi-xt3d": "cannot be used with XT3D",
    "swi-hfb": "cannot be used with the HFB Package",
    "swi2f-grid": "cell bottoms, and cell areas",
    "swi2f-sy": "Specific yield (SY) must be the same",
    "swi2f-tvg": "A TVG file supplies SALTWATER_HEAD",
}

alphaf = 40.0
alphas = 41.0
name = "swi"


def get_ims(sim):
    return flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=200,
        inner_maximum=300,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-9,
        linear_acceleration="bicgstab",
        backtracking_number=20,
        backtracking_tolerance=1.05,
        backtracking_reduction_factor=0.5,
        backtracking_residual_limit=0.002,
    )


def get_oc(gwf, mname):
    return flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{mname}.bud",
        budgetcsv_filerecord=f"{mname}.bud.csv",
        head_filerecord=f"{mname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )


def build_island(ws, case):
    """One-layer strip: recharged land cells draining to seafloor GHB cells."""
    if case == "swi-tr2ss":
        perioddata = [(10.0, 5, 1.0), (1.0, 5, 1.0)]
    else:
        perioddata = [(10.0, 5, 1.0), (10.0, 5, 1.0)]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=2, perioddata=perioddata)
    get_ims(sim)
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    top = np.array(5 * [10.0] + [-20.0, -40.0, -60.0, -80.0, -100.0])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=10, delr=100.0, delc=1.0, top=top, botm=-150.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
    if case == "swi-tr2ss":
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            iconvert=1,
            ss=1.0e-5,
            sy=0.2,
            transient={0: True},
            steady_state={1: True},
        )
    else:
        sto = flopy.mf6.ModflowGwfsto(
            gwf, iconvert=1, ss=1.0e-5, sy=0.2, transient={0: True}
        )
    if case == "swi-tvs":
        flopy.mf6.ModflowUtltvs(
            sto, perioddata={1: [((0, 0, j), "SY", 0.25) for j in range(5)]}
        )
    flopy.mf6.ModflowGwfswi(gwf, zeta_filerecord=f"{name}.zta")
    ghb = [[(0, 0, j), -top[j] / alphaf, 100.0] for j in (5, 6, 7)]
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=ghb)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=np.array(5 * [0.001] + 5 * [0.0]))
    get_oc(gwf, name)
    return sim


def build_strip(ws, case):
    """Confined five-cell strip between two CHD cells."""
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    if case == "swi-tvg-start":
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(10.0, 5, 1.0)])
    else:
        flopy.mf6.ModflowTdis(sim)
    get_ims(sim)
    newton = None if case == "swi-picard" else "newton"
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions=newton
    )
    if case == "swi-tvg-start":
        # interface -40*1 + 41*0.5 = -19.5 lies inside the cells
        top, botm, h1, h2 = -10.0, -100.0, 1.0, 1.0
    else:
        # the interface (below -20) lies below the cells
        top, botm, h1, h2 = 0.0, -10.0, 1.0, 0.5
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=100.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=h1)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0, xt3doptions=(case == "swi-xt3d"))
    if case == "swi-tvg-start":
        flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1.0e-5, sy=0.2, transient={0: True})
        swi = flopy.mf6.ModflowGwfswi(
            gwf, zeta_filerecord=f"{name}.zta", tvg6_filerecord=f"{name}.swi.tvg"
        )
        swi.tvg.initialize(
            auxiliary=["saltwater_head"], aux={0: [np.full((1, 1, 5), 0.5)]}
        )
    else:
        flopy.mf6.ModflowGwfswi(gwf)
    if case == "swi-hfb":
        flopy.mf6.ModflowGwfhfb(
            gwf, stress_period_data=[[(0, 0, 1), (0, 0, 2), 1.0e-3]]
        )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), h1], [(0, 0, 4), h2]])
    get_oc(gwf, name)
    return sim


def build_two_fluid(ws, case):
    """Two-fluid strip whose saltwater model differs from the freshwater one."""
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(10.0, 5, 1.0)])
    ims = get_ims(sim)
    models = []
    for mname in ("freshwater", "saltwater"):
        is_salt = mname == "saltwater"
        gwf = flopy.mf6.ModflowGwf(
            sim, modelname=mname, save_flows=True, newtonoptions="newton"
        )
        botm = -81.0 if (is_salt and case == "swi2f-grid") else -80.0
        flopy.mf6.ModflowGwfdis(
            gwf, nlay=1, nrow=1, ncol=5, delr=100.0, delc=1.0, top=10.0, botm=botm
        )
        flopy.mf6.ModflowGwfic(gwf, strt=0.0 if is_salt else 0.001)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=10.0)
        sy = 0.25 if (is_salt and case == "swi2f-sy") else 0.2
        flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=sy, transient={0: True})
        if case == "swi2f-tvg" and not is_salt:
            swi = flopy.mf6.ModflowGwfswi(gwf, tvg6_filerecord=f"{mname}.swi.tvg")
            swi.tvg.initialize(
                auxiliary=["saltwater_head"], aux={0: [np.full((1, 1, 5), 0.0)]}
            )
        else:
            flopy.mf6.ModflowGwfswi(gwf)
        flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=[[(0, 0, 0), 0.0], [(0, 0, 4), 0.0]]
        )
        get_oc(gwf, mname)
        models.append(gwf)
    flopy.mf6.ModflowSwiswi(
        sim, exgtype="SWI6-SWI6", exgmnamea="freshwater", exgmnameb="saltwater"
    )
    sim.register_ims_package(ims, [m.name for m in models])
    return sim


def build_models(idx, test):
    case = cases[idx]
    if case in ("swi-tvs", "swi-tr2ss"):
        sim = build_island(test.workspace, case)
    elif case.startswith("swi2f"):
        sim = build_two_fluid(test.workspace, case)
    else:
        sim = build_strip(test.workspace, case)
    return sim, None


def check_output(idx, test):
    case = cases[idx]
    ws = test.workspace
    if xfail[idx]:
        lst = (ws / "mfsim.lst").read_text()
        assert errors[case] in lst, f"expected error '{errors[case]}' not found"
        return

    df = pd.read_csv(ws / f"{name}.bud.csv")
    imbalance = (df["TOTAL_IN"] - df["TOTAL_OUT"]).abs().max()
    peak = max(df["TOTAL_IN"].max(), df["TOTAL_OUT"].max())
    assert imbalance < 1.0e-8 + 1.0e-6 * peak, f"budget imbalance {imbalance:.3e}"

    if case == "swi-tr2ss":
        steady = df.iloc[5:]
        swi = steady["STORAGE(SWI)_IN"] + steady["STORAGE(SWI)_OUT"]
        assert np.allclose(swi, 0.0), f"SWI storage in steady state: {swi.values}"
        transient = df.iloc[:5]
        assert (transient["STORAGE(SWI)_IN"] > 0.0).all(), "no transient storage"

    elif case == "swi-tvg-start":
        swi = df["STORAGE(SWI)_IN"] + df["STORAGE(SWI)_OUT"]
        assert np.allclose(swi, 0.0, atol=1.0e-10), (
            f"SWI storage with a consistent starting head: {swi.values}"
        )

    elif case == "swi-picard":
        # confined, fully fresh: Darcy flow K*A*dh/L along four 100-m spans
        qexact = 10.0 * 10.0 * 1.0 * 0.5 / 400.0
        qchd = df["CHD(CHD_0)_IN"].iloc[-1]
        assert abs(qchd - qexact) < 1.0e-6, f"flow {qchd} differs from {qexact}"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
        xfail=xfail[idx],
    )
    test.run()

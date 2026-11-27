"""
Regression test for the MILU0 pivot breakdown in VSF Newton runs.

An infiltration front from a recharge bump moves into a deep, dry Haverkamp
sand (Celia et al., 1990) above an aquifer with fixed heads at the sides.
At the front kr changes by orders of magnitude between neighbouring cells,
and the upstream-kr derivatives make rows of the Newton matrix strongly
non-diagonally dominant. The relaxed modified factorisation (MILU0,
RELAXATION_FACTOR 0.97) then produced pivots of the wrong sign. IMS boosted
the diagonal and finally replaced the pivots with a small value, which left
a useless preconditioner: the inner solves stopped at INNER_MAXIMUM, heads
ran away and the run was killed by a floating-point trap.

IMS now redoes a failed relaxed factorisation without relaxation before it
modifies the diagonal. The test checks that the run terminates normally,
that no pivot was replaced, and that the volume balance closes.

Units are cm and s.
"""

import flopy
import numpy as np
from gwf_test_utils import get_balance_error

DAY = 86400.0
KS, POROSITY, THETA_R = 9.44e-3, 0.287, 0.075
HVK = dict(alphahvk=0.027074, nhvk=3.96, betahvk=0.052408, khvk=4.74)
N, NFINE, NCOARSE = 4, 20, 4
DZFINE, DZCOARSE, TOP, WT_DEPTH = 50.0, 500.0, 10000.0, 900.0


def build_model(ws, exe):
    name = "milupiv"
    dz = np.array([DZFINE] * NFINE + [DZCOARSE] * NCOARSE)
    botm = TOP - np.cumsum(dz)
    zc = TOP - (np.cumsum(dz) - 0.5 * dz)
    nlay = NFINE + NCOARSE
    wt0 = TOP - WT_DEPTH

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="SECONDS",
        nper=2,
        perioddata=[(DAY, 1, 1.0), (4 * DAY, 1, 1.0)],
    )
    tdis.ats.initialize(
        maxats=2,
        perioddata=[(i, 30.0, 1.0, 6 * 3600.0, 2.0, 2.0) for i in range(2)],
        filename=f"{name}.ats",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_dvclose=0.01,
        outer_maximum=100,
        under_relaxation_kappa=0.2,
        backtracking_number=0,
        inner_maximum=300,
        inner_dvclose=0.01,
        rcloserecord=[10.0, "strict"],
        linear_acceleration="BICGSTAB",
        relaxation_factor=0.97,
        preconditioner_levels=0,
        preconditioner_drop_tolerance=0.0,
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=N,
        ncol=N,
        delr=1000.0,
        delc=1000.0,
        top=TOP,
        botm=botm,
    )
    # soil in the unsaturated layers, an anisotropic aquifer below
    k = np.full((nlay, N, N), 0.01)
    k[:NFINE] = KS
    k33 = 0.1 * k
    k33[:NFINE] = KS
    flopy.mf6.ModflowGwfic(gwf, strt=wt0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k, k33=k33)
    flopy.mf6.ModflowGwfsto(gwf, ss=1e-7, sy=0.0, iconvert=0, transient={0: True})
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        kr_averaging="arithmetic",
        porosity=POROSITY,
        satres=THETA_R / POROSITY,
        soil_model="haverkamp",
        **HVK,
    )
    # recharge bump in the middle: 0.2 Ks on average for one day, then
    # 0.001 Ks; the lateral contrast makes the front sharp and uneven
    x = (np.arange(N) + 0.5) / N - 0.5
    bump = 0.4 + 4.0 * np.exp(-(x[:, None] ** 2 + x[None, :] ** 2) / 0.05)
    flopy.mf6.ModflowGwfrcha(gwf, recharge={0: 0.2 * KS * bump, 1: 0.001 * KS * bump})
    edge = [
        (i, j) for i in range(N) for j in range(N) if i in (0, N - 1) or j in (0, N - 1)
    ]
    chd = [[(kk, i, j), wt0] for kk in range(nlay) if zc[kk] < wt0 for (i, j) in edge]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd})
    return sim


def test_milu_pivot_breakdown(function_tmpdir, targets):
    sim = build_model(function_tmpdir, targets["mf6"])
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff)

    lst = (function_tmpdir / "mfsim.lst").read_text()
    assert "ADDED SMALL VALUE TO PIVOT" not in lst, (
        "a relaxed factorisation fell back to replacing pivots"
    )

    err = get_balance_error(function_tmpdir / "milupiv.lst")
    assert abs(err) < 1e-3, f"cumulative volume balance error {err} %"

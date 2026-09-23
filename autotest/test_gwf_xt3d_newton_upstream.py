"""
Test the upstream-saturation treatment of the XT3D Newton formulation.

With the Newton formulation, XT3D scales the flow across a horizontal
connection by the saturation of the upstream cell.  The flow across an XT3D
connection is generally not zero when the heads of the two cells are equal,
because the flow also depends on the heads of the neighboring cells, so the
upstream cell must be selected from the direction of the flow rather than
from the head difference.  Otherwise the product of flow and saturation, and
therefore the discretized residual, jumps when the head ordering of two cells
with different saturations changes.  In a large model this caused a single
cell to flip between two states on every outer iteration and prevented
convergence.

Case 0 (xt3dnr_cont) checks the continuity of the residual directly.  A
small unconfined model with rotated horizontal anisotropy is loaded through
the MODFLOW API, the head of one cell is swept across the head of its
neighbor, which has a different bottom elevation and therefore a different
saturation, and the discretized residual r = A x - b is evaluated at each
head.  The residual of both cells must vary smoothly across the crossing.

Case 1 (xt3dnr_conv) is a small steady-state model in which two rows have
nearly equal transmissivity but very different saturation.  It must converge
with the Newton formulation, its budget must close, and its heads must agree
closely with the standard (non-Newton) XT3D formulation.
"""

import re
from pathlib import Path

import flopy
import numpy as np
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["xt3dnr_cont", "xt3dnr_conv"]

# case 0 grid: the deep cell n and its shallow east neighbor m in the
# middle row are the pair whose head ordering is swept
nlay, nrow, ncol = 1, 3, 5
delr = delc = 100.0
top = 100.0
n_cell = (0, 1, 2)
m_cell = (0, 1, 3)
h_ref = 80.0


def boundary_heads():
    """Heads with gradients in both directions so the XT3D flow across the
    n-m face is not zero when the two cells have the same head."""
    yc = (nrow - 1 - np.arange(nrow)) * delc + delc / 2.0
    xc = np.arange(ncol) * delr + delr / 2.0
    h = np.zeros((nrow, ncol))
    for i in range(nrow):
        for j in range(ncol):
            h[i, j] = h_ref + 0.05 * (yc[i] - yc[1]) + 0.01 * (xc[j] - xc[2])
    return h


def build_cont(ws, exe):
    botm = np.zeros((nlay, nrow, ncol))
    botm[m_cell] = 60.0  # shallow neighbor: saturation 0.5 at h = 80
    sim = flopy.mf6.MFSimulation(sim_name=cases[0], sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_maximum=500,
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        rcloserecord=[1e-8, "STRICT"],
        linear_acceleration="BICGSTAB",
        no_ptcrecord="ALL",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=cases[0], newtonoptions="NEWTON UNDER_RELAXATION"
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
    flopy.mf6.ModflowGwfic(gwf, strt=h_ref)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=1, k=1.0, k22=0.1, angle1=30.0, xt3doptions=True
    )
    hb = boundary_heads()
    spd = [
        [(0, i, j), hb[i, j]]
        for i in range(nrow)
        for j in range(ncol)
        if i in (0, nrow - 1) or j in (0, ncol - 1)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{cases[0]}.hds",
        budget_filerecord=f"{cases[0]}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_conv(ws, exe, newton=True):
    nrow_b, ncol_b = 4, 10
    botm = np.zeros((1, nrow_b, ncol_b))
    k = np.ones((1, nrow_b, ncol_b))
    # rows 1 and 2 have the same transmissivity at h = 90 (1 * 90 and
    # 9 * 10) but saturations of 0.9 and 0.5
    botm[0, 2, :] = 80.0
    k[0, 2, :] = 9.0
    name = cases[1] if newton else cases[1] + "_std"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        outer_maximum=100,
        outer_dvclose=1e-6,
        inner_dvclose=1e-8,
        rcloserecord=[1e-6, "STRICT"],
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="NEWTON UNDER_RELAXATION" if newton else None,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=nrow_b,
        ncol=ncol_b,
        delr=100.0,
        delc=100.0,
        top=100.0,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=90.0)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=1, k=k, k22=0.1, angle1=30.0, xt3doptions=True
    )
    spd = [[(0, i, 0), 92.0] for i in range(nrow_b)] + [
        [(0, i, ncol_b - 1), 88.0] for i in range(nrow_b)
    ]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    if idx == 0:
        sim = build_cont(test.workspace, "libmf6")
    else:
        sim = build_conv(test.workspace, "mf6")
    return sim, None


def residual_sweep(libmf6, ws):
    """Sweep the head of cell n across the head of cell m and return the
    discretized residuals of both cells."""
    from modflowapi import ModflowApi

    mf6 = ModflowApi(libmf6, working_directory=str(ws))
    mf6.initialize()
    mf6.prepare_time_step(mf6.get_time_step())
    mf6.prepare_solve()
    xtag = mf6.get_var_address("X", "SLN_1")
    atag = mf6.get_var_address("AMAT", "SLN_1")
    btag = mf6.get_var_address("RHS", "SLN_1")
    ia = mf6.get_value(mf6.get_var_address("IA", "SLN_1")) - 1
    ja = mf6.get_value(mf6.get_var_address("JA", "SLN_1")) - 1
    nn = n_cell[1] * ncol + n_cell[2]
    mm = m_cell[1] * ncol + m_cell[2]
    base = boundary_heads().ravel()
    base[mm] = h_ref
    deltas = np.linspace(-0.05, 0.05, 41)
    res = np.zeros((deltas.size, 2))
    for k, d in enumerate(deltas):
        x = base.copy()
        x[nn] = h_ref + d
        mf6.set_value(xtag, x.copy())
        # solve() assembles A and b for the heads just set (the Newton
        # terms cancel in A x - b) and then overwrites x with the linear
        # solution, which is discarded
        mf6.solve()
        a = mf6.get_value(atag)
        b = mf6.get_value(btag)
        for irow, icell in ((0, nn), (1, mm)):
            cols = ja[ia[icell] : ia[icell + 1]]
            res[k, irow] = a[ia[icell] : ia[icell + 1]] @ x[cols] - b[icell]
    mf6.finalize_solve()
    mf6.finalize_time_step()
    mf6.finalize()
    return deltas, res


def api_func(exe, idx, model_ws):
    """Check that the residual of both cells is continuous across the point
    where the heads of the two cells are equal."""
    deltas, res = residual_sweep(exe, model_ws)
    ok = True
    buff = []
    for icell, label in enumerate(("n", "m")):
        dr = np.abs(np.diff(res[:, icell]))
        ratio = dr.max() / np.median(dr)
        buff.append(
            f"cell {label}: median step {np.median(dr):.4g}, "
            f"largest step {dr.max():.4g}, ratio {ratio:.2f}"
        )
        # the residual is nearly linear in the swept head, so consecutive
        # steps should be nearly equal; a discontinuity at the crossing
        # shows up as one step that is orders of magnitude larger
        if ratio > 2.0:
            ok = False
            buff.append(f"  residual of cell {label} is discontinuous")
    print("\n".join(buff))
    return ok, buff


def check_output(idx, test):
    ws = Path(test.workspace)
    if idx == 0:
        return
    name = cases[1]
    lst = (ws / "mfsim.lst").read_text()
    nouter = len(
        re.findall(r"^\s+(?:Model|Under-relaxation|Backtracking)\s+\d+", lst, re.M)
    )
    print(f"Newton XT3D outer iterations: {nouter}")
    assert nouter < 60, f"Newton XT3D needed {nouter} outer iterations"
    disc = [
        abs(float(v))
        for v in re.findall(
            r"PERCENT DISCREPANCY =\s+(-?\d+\.\d+)", (ws / f"{name}.lst").read_text()
        )
    ]
    assert disc and max(disc) < 0.01, f"budget discrepancy {disc}"
    # reference: same model with the standard formulation
    ws_std = ws / "std"
    sim = build_conv(ws_std, test.targets["mf6"], newton=False)
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert success, "standard-formulation reference run failed"
    h_newton = flopy.utils.HeadFile(ws / f"{name}.hds").get_data()
    h_std = flopy.utils.HeadFile(ws_std / f"{name}_std.hds").get_data()
    dmax = np.abs(h_newton - h_std).max()
    print(f"max |Newton - standard| head difference: {dmax:.4f}")
    # the two formulations discretize saturated thickness slightly
    # differently, so they agree closely but not exactly
    assert dmax < 0.2, f"Newton and standard XT3D heads differ by {dmax}"


@requires_pkg("modflowapi")
def test_residual_continuity(function_tmpdir, targets):
    idx = 0
    test = TestFramework(
        name=cases[idx],
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()


def test_convergence(function_tmpdir, targets):
    idx = 1
    test = TestFramework(
        name=cases[idx],
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

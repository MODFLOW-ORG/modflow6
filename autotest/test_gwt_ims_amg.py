"""
Integration test for the AMG preconditioner applied to a GWT (groundwater
transport) model with upstream-weighted advection-dispersion.

The GWT matrix is strongly nonsymmetric because upstream weighting adds an
advective flux to a neighbour's matrix entry only when flow comes from that
neighbour, so a_ij != a_ji.  The difference grows with the cell Peclet number
Pe = v*dx / (2*D), so a higher Pe makes the matrix more nonsymmetric and
stresses AMG's handling of nonsymmetric systems.

Model setup:
  - Grid: nlay=2, nrow=15, ncol=25 (750 active cells), dx=dy=100 m, dz=10 m
  - Steady GWF: west CHD=10 m, east CHD=0 m; K=[50, 25] m/d per layer with a
    checkerboard multiplier (x 2 or x 0.5) to add lateral heterogeneity
  - GWF IMS: CG + adaptive omega (RELAX=0) + ILU (5 levels) — symmetric GWF
  - Transient GWT: nper=3, perlen=1000 d, nstp=5; upstream advection;
    CNC west face C=1.0; initial C=0; porosity=0.25
  - GWT IMS: BiCGSTAB + fixed omega (RELAX=1); AMG or ILU preconditioner

Two longitudinal dispersivity cases:
  - gwt_amg_lope  (alpha_L=20 m, Pe = v·dx/(2·alpha_L·v) = dx/(2·alpha_L)
                    = 100/(40) = 2.5): moderate nonsymmetry
  - gwt_amg_hipe  (alpha_L= 6 m, Pe = 100/12 ≈ 8.3): strong nonsymmetry

Each case runs an AMG sim (GWT solved with AMG) and an ILU reference sim
(GWT solved with MILUT).  Final concentrations are compared within ctol=1e-4.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from amg_test_utils import inner_iters_per_timestep
from framework import TestFramework

paktest = "ims"
cmp_prefix = "mf6"

cases = [
    "gwt_amg_lope",  # alpha_L=20 m, Pe≈2.5
    "gwt_amg_hipe",  # alpha_L= 6 m, Pe≈8.3
]
alpha_L_values = [20.0, 6.0]

# ---------------------------------------------------------------------------
# Grid and flow model parameters
# ---------------------------------------------------------------------------

nlay, nrow, ncol = 2, 15, 25
delr, delc = 100.0, 100.0  # m
top = 20.0
botm = [10.0, 0.0]

chd_west = 10.0
chd_east = 0.0
k_base = np.array([50.0, 25.0])  # m/d per layer

# Transport
nper = 3
perlen = 1000.0  # d
nstp = 5
cnc_conc = 1.0
init_conc = 0.0
porosity = 0.25

# Solver levels
amg_levels = 5
ilu_levels = 5

ctol = 1.0e-4  # concentration tolerance for AMG-vs-ILU comparison


# ---------------------------------------------------------------------------
# K field
# ---------------------------------------------------------------------------


def _build_k_field():
    """Two-layer K array with a checkerboard multiplier (x 2 or x 0.5)."""
    k = np.empty((nlay, nrow, ncol))
    for ilay in range(nlay):
        for irow in range(nrow):
            for icol in range(ncol):
                mult = 2.0 if (irow + icol) % 2 == 0 else 0.5
                k[ilay, irow, icol] = k_base[ilay] * mult
    return k


# ---------------------------------------------------------------------------
# Model-building helpers
# ---------------------------------------------------------------------------


def _build_gwf(sim, gwf_name):
    """Steady GWF model with west/east CHD boundaries."""
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{gwf_name}.dis",
    )

    strt = np.broadcast_to(
        np.linspace(chd_west, chd_east, ncol), (nlay, nrow, ncol)
    ).copy()
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwf_name}.ic")

    k = _build_k_field()
    flopy.mf6.ModflowGwfnpf(
        gwf, k=k, save_specific_discharge=True, filename=f"{gwf_name}.npf"
    )

    spd = [[(ilay, irow, 0), chd_west] for ilay in range(nlay) for irow in range(nrow)]
    spd += [
        [(ilay, irow, ncol - 1), chd_east]
        for ilay in range(nlay)
        for irow in range(nrow)
    ]
    flopy.mf6.modflow.ModflowGwfchd(gwf, stress_period_data=spd)

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwf_name}.hds",
        budget_filerecord=f"{gwf_name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return gwf


def _build_gwt(sim, gwf_name, gwt_name, alpha_L):
    """Transient GWT model with upstream-weighted advection-dispersion."""
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwt_name,
        model_nam_file=f"{gwt_name}.nam",
    )
    gwt.name_file.save_flows = True

    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{gwt_name}.dis",
    )
    flopy.mf6.ModflowGwtic(gwt, strt=init_conc, filename=f"{gwt_name}.ic")

    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwt_name}.adv")
    flopy.mf6.ModflowGwtdsp(
        gwt,
        alh=alpha_L,
        ath1=alpha_L * 0.1,
        diffc=0.0,
        filename=f"{gwt_name}.dsp",
    )
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwt_name}.mst")

    # SSM required when GWF has boundary packages; no auxiliary sources here
    # because CNC directly fixes concentrations at the west boundary cells.
    flopy.mf6.ModflowGwtssm(gwt, filename=f"{gwt_name}.ssm")

    # Constant concentration on west face (all layers)
    cnc_spd = [
        [(ilay, irow, 0), cnc_conc] for ilay in range(nlay) for irow in range(nrow)
    ]
    flopy.mf6.ModflowGwtcnc(gwt, stress_period_data={0: cnc_spd})

    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwt_name}.ucn",
        saverecord=[("CONCENTRATION", "ALL")],
    )

    # Couple GWF flow field to GWT
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_name,
        exgmnameb=gwt_name,
        filename=f"{gwt_name}.gwfgwt",
    )
    return gwt


def _gwf_ims_kwargs(gwf_name):
    return dict(
        print_option="SUMMARY",
        outer_dvclose=1.0e-8,
        outer_maximum=100,
        inner_dvclose=1.0e-12,
        inner_maximum=300,
        rcloserecord="1.0e-8 strict",
        linear_acceleration="cg",
        preconditioner_levels=ilu_levels,
        relaxation_factor=0.0,
        filename=f"{gwf_name}_ims.ims",
    )


def _gwt_ims_kwargs(gwt_name, levels):
    return dict(
        print_option="SUMMARY",
        outer_dvclose=1.0e-8,
        outer_maximum=100,
        inner_dvclose=1.0e-12,
        inner_maximum=300,
        rcloserecord="1.0e-8 strict",
        linear_acceleration="bicgstab",
        preconditioner_levels=levels,
        relaxation_factor=1.0,
        filename=f"{gwt_name}_ims.ims",
        csv_outer_output_filerecord=f"{gwt_name}_outer.csv",
    )


def _build_sim(name, ws, alpha_L, use_amg):
    """Build a GWF+GWT simulation; use_amg enables AMG on the GWT IMS."""
    gwf_name = f"{name}_gwf"
    gwt_name = name

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        nper=nper,
        perioddata=[(perlen, nstp, 1.0)] * nper,
    )

    gwf_ims = flopy.mf6.ModflowIms(sim, **_gwf_ims_kwargs(gwf_name))
    gwt_ims = flopy.mf6.ModflowIms(
        sim,
        **_gwt_ims_kwargs(gwt_name, amg_levels if use_amg else ilu_levels),
        preconditioner_type="AMG" if use_amg else None,
    )

    gwf = _build_gwf(sim, gwf_name)
    gwt = _build_gwt(sim, gwf_name, gwt_name, alpha_L)

    sim.register_ims_package(gwf_ims, [gwf.name])
    sim.register_ims_package(gwt_ims, [gwt.name])

    sim.write_simulation()
    return sim


# ---------------------------------------------------------------------------
# TestFramework interface
# ---------------------------------------------------------------------------


def build_models(idx, test):
    name = cases[idx]
    alpha_L = alpha_L_values[idx]
    cmp_ws = Path(test.workspace) / cmp_prefix

    amg_sim = _build_sim(name, test.workspace, alpha_L, use_amg=True)
    ilu_sim = _build_sim(name, cmp_ws, alpha_L, use_amg=False)
    return amg_sim, ilu_sim


def check_output(idx, test):
    name = cases[idx]

    conc_amg = flopy.utils.HeadFile(
        Path(test.workspace) / f"{name}.ucn",
        precision="double",
        text="CONCENTRATION",
    ).get_alldata()

    conc_ilu = flopy.utils.HeadFile(
        Path(test.workspace) / cmp_prefix / f"{name}.ucn",
        precision="double",
        text="CONCENTRATION",
    ).get_alldata()

    assert conc_amg.shape == conc_ilu.shape, (
        f"Concentration array shapes differ: AMG {conc_amg.shape} vs "
        f"ILU {conc_ilu.shape}"
    )

    max_diff = np.max(np.abs(conc_amg - conc_ilu))
    assert max_diff < ctol, (
        f"AMG and ILU concentrations differ by {max_diff:.3e} (tolerance {ctol:.3e})"
    )

    amg_ts = inner_iters_per_timestep(Path(test.workspace) / f"{name}_outer.csv")
    ilu_ts = inner_iters_per_timestep(
        Path(test.workspace) / cmp_prefix / f"{name}_outer.csv"
    )
    assert amg_ts.sum() < ilu_ts.sum(), (
        f"AMG total inner iterations ({amg_ts.sum()}) >= ILU ({ilu_ts.sum()}) "
        f"summed over {len(amg_ts)} time steps"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        overwrite=False,
        compare=None,
    )
    test.run()

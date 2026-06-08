"""
Tests the AUXLENGTHNAME option for the WEL package AUTO_FLOW_REDUCE feature.

Two adjacent unconfined cells (10 m thick, hk=1e-4 m/d so effectively isolated)
are each pumped at -0.25 m³/d with FLOW_REDUCTION_LENGTH and AUXLENGTHNAME active.
The per-well LENGTH auxiliary variable assigns different thresholds:
  Cell 1 (col 0): LENGTH = 2.0 m  →  tp = 0 + 2.0 = 2.0 m
  Cell 2 (col 1): LENGTH = 5.0 m  →  tp = 0 + 5.0 = 5.0 m

Cell 2 starts reducing pumping earlier (when h ≤ 5 m, at t ≈ 2 d), retaining more
water, so h_cell2 > h_cell1 at t = 4 d.

Main model: Newton + AUXLENGTHNAME (per-well threshold via aux variable)
Reference  (mf6/): Newton + two separate WEL packages (one per cell), each with
  an explicit global AUTO_FLOW_REDUCE value equal to the per-well aux value.

The TestFramework compare="mf6" verifies that the AUXLENGTHNAME model gives
identical results (within 1 mm) to the explicit two-package reference, proving
that AUXLENGTHNAME correctly overrides the global threshold per well.
check_output() additionally checks the physics: h_cell2 > h_cell1.
"""

import pathlib
import subprocess

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["wel03"]

# --- Grid: two adjacent effectively-isolated cells, 10 m thick ---
nlay, nrow, ncol = 1, 1, 2
top = 10.0
botm = [0.0]
delr = [1.0, 1.0]
delc = 1.0

# --- Hydraulic properties ---
hk = 1e-4  # very low: isolates the two cells
sy = 0.1
strt = 10.0

# --- WEL parameters ---
wellq = -0.25  # m³/d
# Per-well FLOW_REDUCTION_LENGTH values (via AUXLENGTHNAME aux variable):
_len_cell0 = 2.0  # tp = botm + 2.0 = 2.0 m for cell 0
_len_cell1 = 5.0  # tp = botm + 5.0 = 5.0 m for cell 1

# --- Time: 4 days, 40 equal steps of 0.1 d ---
nper = 1
tdis_rc = [(4.0, 40, 1.0)]

# --- Solver ---
nouter, ninner = 500, 300
hclose = 1e-9
rcloserecord = 1e-3
relax = 1.0


def _build_gwf_base(sim, name):
    """Create the GWF model with DIS, IC, NPF, STO, OC; no WEL."""
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="NEWTON UNDER_RELAXATION",
        save_flows=True,
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
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk, k33=hk)
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=sy,
        transient={0: True},
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return gwf


def _build_tdis_ims(sim, name):
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rcloserecord,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )


def _add_auxlengthname(wel_file, aux_name):
    """Inject AUXLENGTHNAME keyword into an already-written .wel OPTIONS block."""
    content = wel_file.read_text()
    content = content.replace(
        "END options",
        f"  AUXLENGTHNAME {aux_name}\nEND options",
        1,
    )
    wel_file.write_text(content)


def build_models(idx, test):
    name = cases[idx]

    # --- Main model: Newton + AUXLENGTHNAME ---
    # Flopy does not yet know about AUXLENGTHNAME, so we write the model and
    # then inject the new keyword into the OPTIONS block of the .wel file.
    sim_main = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    _build_tdis_ims(sim_main, name)
    gwf_main = _build_gwf_base(sim_main, name)

    # WEL with AUXILIARY LENGTH and AUTO_FLOW_REDUCE + FLOW_REDUCTION_LENGTH.
    # The global AUTO_FLOW_REDUCE value is a placeholder (1.0); it is overridden
    # per-well by AUXLENGTHNAME once the keyword is injected below.
    flopy.mf6.ModflowGwfwel(
        gwf_main,
        auxiliary=["LENGTH"],
        auto_flow_reduce=1.0,
        flow_reduction_length=True,
        print_flows=True,
        stress_period_data={
            0: [
                [(0, 0, 0), wellq, _len_cell0],
                [(0, 0, 1), wellq, _len_cell1],
            ]
        },
    )

    sim_main.write_simulation(silent=True)
    # Inject AUXLENGTHNAME into the written .wel file
    _add_auxlengthname(pathlib.Path(test.workspace) / f"{name}.wel", "LENGTH")

    # Run the main model now; the framework never runs it because we return None.
    result = subprocess.run(
        [str(test.targets["mf6"])],
        cwd=test.workspace,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Main AUXLENGTHNAME model failed:\n{result.stdout}\n{result.stderr}"
    )

    # --- Reference model (mf6/): Newton + two explicit WEL packages ---
    # Each package covers one cell and carries its own AUTO_FLOW_REDUCE value,
    # so the results should be bit-for-bit identical to the AUXLENGTHNAME model.
    ref_ws = pathlib.Path(test.workspace) / "mf6"
    sim_ref = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=str(ref_ws)
    )
    _build_tdis_ims(sim_ref, name)
    gwf_ref = _build_gwf_base(sim_ref, name)

    flopy.mf6.ModflowGwfwel(
        gwf_ref,
        pname="WEL-1",
        filename="wel1.wel",
        auto_flow_reduce=_len_cell0,
        flow_reduction_length=True,
        print_flows=True,
        stress_period_data={0: [[(0, 0, 0), wellq]]},
    )
    flopy.mf6.ModflowGwfwel(
        gwf_ref,
        pname="WEL-2",
        filename="wel2.wel",
        auto_flow_reduce=_len_cell1,
        flow_reduction_length=True,
        print_flows=True,
        stress_period_data={0: [[(0, 0, 1), wellq]]},
    )

    # Return None for main (already written); framework writes sim_ref.
    return None, sim_ref


def check_output(idx, test):
    name = cases[idx]

    h_main = flopy.utils.HeadFile(
        pathlib.Path(test.workspace) / f"{name}.hds"
    ).get_alldata()
    h_ref = flopy.utils.HeadFile(
        pathlib.Path(test.workspace) / "mf6" / f"{name}.hds"
    ).get_alldata()

    # AUXLENGTHNAME model and explicit two-package reference must agree within 1 mm.
    max_diff = float(np.max(np.abs(h_main - h_ref)))
    assert max_diff < 1e-3, (
        f"{name}: max head difference between AUXLENGTHNAME and "
        f"two-package reference = {max_diff:.2e} m (must be < 1e-3 m)"
    )

    # Physics check: cell with tp=5.0 m (col 1) reduces pumping earlier and
    # retains more water than cell with tp=2.0 m (col 0) at t = 4 d.
    h_final = h_main[-1, 0, 0, :]  # shape (ncol,)
    h_cell0 = float(h_final[0])  # tp = 2.0 m
    h_cell1 = float(h_final[1])  # tp = 5.0 m
    assert h_cell1 > h_cell0, (
        f"{name}: cell with tp=5.0 m should have higher head than cell with "
        f"tp=2.0 m, but got h_cell0={h_cell0:.4f} m, h_cell1={h_cell1:.4f} m"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare="mf6",
    )
    test.run()

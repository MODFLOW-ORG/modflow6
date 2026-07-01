"""
Tests the AUTO_FLOW_REDUCE_AUXNAME option for the WEL package AUTO_FLOW_REDUCE feature.

Two unconfined cells (10 m thick) separated by an inactive cell -- so they are
truly independent 1-cell problems -- are each pumped at -0.25 m³/d with
FLOW_REDUCTION_LENGTH and AUTO_FLOW_REDUCE_AUXNAME active. The per-well LENGTH auxiliary
variable assigns different thresholds:
  Cell 0 (col 0): LENGTH = 2.0 m  →  tp = 0 + 2.0 = 2.0 m
  Cell 1 (col 2): LENGTH = 5.0 m  →  tp = 0 + 5.0 = 5.0 m

Cell 1 starts reducing pumping earlier (when h ≤ 5 m, at t ≈ 2 d), retaining more
water, so h_cell1 > h_cell0 at t = 4 d.

Main model: Newton + AUTO_FLOW_REDUCE_AUXNAME (per-well threshold via aux variable)
Reference  (mf6/): Newton + two separate WEL packages (one per cell), each with
  an explicit global AUTO_FLOW_REDUCE value equal to the per-well aux value.
Standard   (std/): standard formulation + AUTO_FLOW_REDUCE_AUXNAME (otherwise
  identical to the main model).

The TestFramework compare="mf6" verifies that the AUTO_FLOW_REDUCE_AUXNAME model gives
identical results (within 1 mm) to the explicit two-package reference, proving
that AUTO_FLOW_REDUCE_AUXNAME correctly overrides the global threshold per well.
check_output() additionally:
  - checks that the standard and Newton formulations give identical results
    (the AUTO_FLOW_REDUCE reduction is applied in wel_cf for both); and
  - checks the physics: h_cell1 > h_cell0.

test_afr_auxname_messages() separately exercises the option's warning and
error paths (AUTO_FLOW_REDUCE_AUXNAME without AUTO_FLOW_REDUCE or FLOW_REDUCTION_LENGTH,
a missing aux name, and no aux variables at all).
"""

import pathlib
import subprocess

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["wel03"]

# --- Grid: two pumped cells separated by an inactive cell, 10 m thick ---
# The two pumped cells (col 0 and col 2) share no cell face, so they are truly
# independent 1-cell problems. This structural isolation (rather than a very
# low K between adjacent cells) is important: with only a near-zero coupling the
# iterative solver can transiently overshoot one cell below its bottom, and the
# standard formulation then irreversibly converts it to dry (HDRY) -- which
# cell depends on the linear accelerator. With a fully inactive separator the
# Newton and standard formulations converge to identical heads.
nlay, nrow, ncol = 1, 1, 3
top = 10.0
botm = [0.0]
delr = [1.0, 1.0, 1.0]
delc = 1.0
idomain = [[[1, 0, 1]]]  # middle cell inactive

# --- Cell indices of the two pumped cells ---
_col_cell0 = 0
_col_cell1 = 2

# --- Hydraulic properties ---
hk = 1.0  # isolation is structural (inactive separator), so K is immaterial
sy = 0.1
strt = 10.0

# --- WEL parameters ---
wellq = -0.25  # m³/d
# Per-well FLOW_REDUCTION_LENGTH values (via AUTO_FLOW_REDUCE_AUXNAME aux variable):
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


def _build_gwf_base(sim, name, newton=True):
    """Create the GWF model with DIS, IC, NPF, STO, OC; no WEL.

    When newton is True the model uses the Newton-Raphson formulation;
    when False it uses the standard formulation. The AUTO_FLOW_REDUCE
    pumping reduction is applied in wel_cf regardless of formulation, so
    both should converge to the same heads.
    """
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions="NEWTON UNDER_RELAXATION" if newton else None,
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
        idomain=idomain,
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


def _add_afr_auxname(wel_file, aux_name):
    """Inject AUTO_FLOW_REDUCE_AUXNAME into an already-written .wel OPTIONS block."""
    content = wel_file.read_text()
    content = content.replace(
        "END options",
        f"  AUTO_FLOW_REDUCE_AUXNAME {aux_name}\nEND options",
        1,
    )
    wel_file.write_text(content)


def _run_mf6(exe, ws):
    """Run the mf6 executable in ws and return the CompletedProcess."""
    return subprocess.run([str(exe)], cwd=str(ws), capture_output=True, text=True)


def _build_auxlength_sim(ws, name, newton):
    """Build the per-well AUTO_FLOW_REDUCE_AUXNAME model (DIS/IC/NPF/STO/OC + WEL).

    The global AUTO_FLOW_REDUCE value is a placeholder (1.0); it is overridden
    per-well by AUTO_FLOW_REDUCE_AUXNAME once the keyword is injected into the
    .wel file.
    """
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=str(ws)
    )
    _build_tdis_ims(sim, name)
    gwf = _build_gwf_base(sim, name, newton=newton)
    flopy.mf6.ModflowGwfwel(
        gwf,
        auxiliary=["LENGTH"],
        auto_flow_reduce=1.0,
        flow_reduction_length=True,
        print_flows=True,
        stress_period_data={
            0: [
                [(0, 0, _col_cell0), wellq, _len_cell0],
                [(0, 0, _col_cell1), wellq, _len_cell1],
            ]
        },
    )
    return sim


def build_models(idx, test):
    name = cases[idx]

    # --- Main model: Newton + AUTO_FLOW_REDUCE_AUXNAME ---
    # Flopy does not yet know about AUTO_FLOW_REDUCE_AUXNAME, so we write the model and
    # then inject the new keyword into the OPTIONS block of the .wel file.
    sim_main = _build_auxlength_sim(test.workspace, name, newton=True)
    sim_main.write_simulation(silent=True)
    # Inject AUTO_FLOW_REDUCE_AUXNAME into the written .wel file
    _add_afr_auxname(pathlib.Path(test.workspace) / f"{name}.wel", "LENGTH")

    # Run the main model now; the framework never runs it because we return None.
    result = _run_mf6(test.targets["mf6"], test.workspace)
    assert result.returncode == 0, (
        f"Main AUTO_FLOW_REDUCE_AUXNAME model failed:\n{result.stdout}\n{result.stderr}"
    )

    # --- Standard-formulation AUTO_FLOW_REDUCE_AUXNAME model (std/) ---
    # Identical to the main model but without the Newton-Raphson formulation.
    # The pumping reduction is applied in wel_cf regardless of formulation, so
    # this must converge to the same heads as the Newton main model.
    std_ws = pathlib.Path(test.workspace) / "std"
    sim_std = _build_auxlength_sim(std_ws, name, newton=False)
    sim_std.write_simulation(silent=True)
    _add_afr_auxname(std_ws / f"{name}.wel", "LENGTH")
    result = _run_mf6(test.targets["mf6"], std_ws)
    assert result.returncode == 0, (
        f"Standard-formulation AUTO_FLOW_REDUCE_AUXNAME model failed:"
        f"\n{result.stdout}\n{result.stderr}"
    )

    # --- Reference model (mf6/): Newton + two explicit WEL packages ---
    # Each package covers one cell and carries its own AUTO_FLOW_REDUCE value,
    # so the results should be bit-for-bit identical to the per-well model.
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
        stress_period_data={0: [[(0, 0, _col_cell0), wellq]]},
    )
    flopy.mf6.ModflowGwfwel(
        gwf_ref,
        pname="WEL-2",
        filename="wel2.wel",
        auto_flow_reduce=_len_cell1,
        flow_reduction_length=True,
        print_flows=True,
        stress_period_data={0: [[(0, 0, _col_cell1), wellq]]},
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
    h_std = flopy.utils.HeadFile(
        pathlib.Path(test.workspace) / "std" / f"{name}.hds"
    ).get_alldata()

    # Compare only active cells (the inactive separator carries 1e30).
    active = h_main < 1e29

    # Per-well model and explicit two-package reference must agree within 1 mm.
    max_diff = float(np.max(np.abs(h_main[active] - h_ref[active])))
    assert max_diff < 1e-3, (
        f"{name}: max head difference between AUTO_FLOW_REDUCE_AUXNAME and "
        f"two-package reference = {max_diff:.2e} m (must be < 1e-3 m)"
    )

    # Standard and Newton formulations must give identical results (within 1 mm):
    # the AUTO_FLOW_REDUCE pumping reduction is applied in wel_cf for both.
    max_diff_form = float(np.max(np.abs(h_main[active] - h_std[active])))
    assert max_diff_form < 1e-3, (
        f"{name}: max head difference between Newton and standard "
        f"formulation = {max_diff_form:.2e} m (must be < 1e-3 m)"
    )

    # Physics check: cell with tp=5.0 m reduces pumping earlier and retains
    # more water than cell with tp=2.0 m at t = 4 d.
    h_final = h_main[-1, 0, 0, :]  # shape (ncol,)
    h_cell0 = float(h_final[_col_cell0])  # tp = 2.0 m
    h_cell1 = float(h_final[_col_cell1])  # tp = 5.0 m
    assert h_cell1 > h_cell0, (
        f"{name}: cell with tp=5.0 m (col {_col_cell1}) should have higher head "
        f"than cell with tp=2.0 m (col {_col_cell0}), but got "
        f"h_cell0={h_cell0:.4f} m, h_cell1={h_cell1:.4f} m"
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


def _run_message_case(exe, ws, afr, frl, aux, inject):
    """Build a minimal one-cell WEL model with the given options, inject an
    AUTO_FLOW_REDUCE_AUXNAME keyword, run it, and return (rc, normalized stdout).

    Parameters control whether AUTO_FLOW_REDUCE (afr) and FLOW_REDUCTION_LENGTH
    (frl) are written and whether an AUX variable (aux) is declared. The
    AUTO_FLOW_REDUCE_AUXNAME value injected into the OPTIONS block is given by inject.
    """
    ws = pathlib.Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    sim = flopy.mf6.MFSimulation(
        sim_name="m", version="mf6", exe_name="mf6", sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=1, delr=1.0, delc=1.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=0.1, transient={0: True})

    kw = {}
    if afr is not None:
        kw["auto_flow_reduce"] = afr
    if frl:
        kw["flow_reduction_length"] = True
    if aux:
        kw["auxiliary"] = aux
        spd = [[(0, 0, 0), -0.1, 1.0]]
    else:
        spd = [[(0, 0, 0), -0.1]]
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: spd}, **kw)

    sim.write_simulation(silent=True)
    _add_afr_auxname(ws / "m.wel", inject)
    result = _run_mf6(exe, ws)
    # Collapse whitespace so wrapped warning/error lines can be matched.
    return result.returncode, " ".join(result.stdout.split())


# (afr, frl, aux, inject, expect_success, required substrings)
_MESSAGE_CASES = [
    pytest.param(
        None,
        False,
        ["LENGTH"],
        "LENGTH",
        True,
        [
            "AUTO_FLOW_REDUCE_AUXNAME is specified but AUTO_FLOW_REDUCE is "
            "not specified",
            "AUTO_FLOW_REDUCE_AUXNAME is specified but FLOW_REDUCTION_LENGTH "
            "is not specified",
        ],
        id="no_afr_no_frl_warns",
    ),
    pytest.param(
        1.0,
        False,
        ["LENGTH"],
        "LENGTH",
        True,
        [
            "AUTO_FLOW_REDUCE_AUXNAME is specified but FLOW_REDUCTION_LENGTH "
            "is not specified",
            "interpreted as a fraction of the cell thickness",
        ],
        id="no_frl_fraction_warns",
    ),
    pytest.param(
        1.0,
        True,
        None,
        "LENGTH",
        False,
        [
            "AUTO_FLOW_REDUCE_AUXNAME was specified as LENGTH but no AUX "
            "variables specified"
        ],
        id="no_aux_errors",
    ),
    pytest.param(
        1.0,
        True,
        ["LENGTH"],
        "BADNAME",
        False,
        [
            "AUTO_FLOW_REDUCE_AUXNAME was specified as BADNAME but no AUX "
            "variable found with this name"
        ],
        id="bad_aux_name_errors",
    ),
]


@pytest.mark.parametrize(
    "afr, frl, aux, inject, expect_success, required", _MESSAGE_CASES
)
def test_afr_auxname_messages(
    afr, frl, aux, inject, expect_success, required, function_tmpdir, targets
):
    """Exercise the AUTO_FLOW_REDUCE_AUXNAME warning and error paths."""
    returncode, out = _run_message_case(
        targets["mf6"], function_tmpdir, afr, frl, aux, inject
    )
    if expect_success:
        assert returncode == 0, f"expected success but mf6 failed:\n{out}"
    else:
        assert returncode != 0, f"expected failure but mf6 succeeded:\n{out}"
    for substring in required:
        assert substring in out, (
            f"expected message not found:\n  {substring!r}\nin output:\n{out}"
        )


# Cell thickness used by the bounds-check model below (top - botm).
_BOUNDS_THICK = 10.0


def _run_bounds_case(exe, ws, frl, aux_values):
    """Build a one-cell WEL model whose AUTO_FLOW_REDUCE_AUXNAME aux variable
    takes the given per-period values, run it, and return (rc, normalized
    stdout).

    aux_values is a list with one value per stress period, so a multi-period
    list exercises the per-period bounds check (the value may change between
    periods). The cell is _BOUNDS_THICK thick.
    """
    ws = pathlib.Path(ws)
    ws.mkdir(parents=True, exist_ok=True)
    nper = len(aux_values)
    sim = flopy.mf6.MFSimulation(
        sim_name="m", version="mf6", exe_name="mf6", sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(1.0, 1, 1.0)] * nper)
    flopy.mf6.ModflowIms(sim)
    gwf = flopy.mf6.ModflowGwf(sim, modelname="m", save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=1,
        delr=1.0,
        delc=1.0,
        top=_BOUNDS_THICK,
        botm=[0.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=_BOUNDS_THICK)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=0.0, sy=0.1, transient={0: True})

    kw = {"auxiliary": ["LENGTH"], "auto_flow_reduce": 1.0}
    if frl:
        kw["flow_reduction_length"] = True
    spd = {p: [[(0, 0, 0), -0.1, av]] for p, av in enumerate(aux_values)}
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data=spd, **kw)

    sim.write_simulation(silent=True)
    _add_afr_auxname(ws / "m.wel", "LENGTH")
    result = _run_mf6(exe, ws)
    return result.returncode, " ".join(result.stdout.split())


_FRAC_ERR = (
    "must be greater than 0 and less than or equal to 1 when "
    "FLOW_REDUCTION_LENGTH is not specified"
)
_LEN_ERR = "must be greater than 0 and less than or equal to the cell thickness"

# (frl, aux_values, expect_success, required substrings)
_BOUNDS_CASES = [
    # Fraction mode: valid range is (0, 1].
    pytest.param(False, [-0.5], False, [_FRAC_ERR], id="frac_negative"),
    pytest.param(False, [0.0], False, [_FRAC_ERR], id="frac_zero"),
    pytest.param(False, [1.5], False, [_FRAC_ERR], id="frac_above_one"),
    pytest.param(False, [0.5], True, [], id="frac_valid"),
    # Length mode: valid range is (0, cell thickness].
    pytest.param(True, [-1.0], False, [_LEN_ERR], id="length_negative"),
    pytest.param(True, [0.0], False, [_LEN_ERR], id="length_zero"),
    pytest.param(
        True,
        [_BOUNDS_THICK + 5.0],
        False,
        [_LEN_ERR],
        id="length_above_thick",
    ),
    pytest.param(True, [5.0], True, [], id="length_valid"),
    # A value valid in period 1 but invalid in period 2 must still be caught,
    # confirming the check runs every stress period.
    pytest.param(
        True,
        [5.0, _BOUNDS_THICK + 5.0],
        False,
        [_LEN_ERR],
        id="transient_becomes_invalid",
    ),
]


@pytest.mark.parametrize("frl, aux_values, expect_success, required", _BOUNDS_CASES)
def test_afr_auxname_bounds(
    frl, aux_values, expect_success, required, function_tmpdir, targets
):
    """Per-well AUTO_FLOW_REDUCE_AUXNAME values out of range must error."""
    returncode, out = _run_bounds_case(targets["mf6"], function_tmpdir, frl, aux_values)
    if expect_success:
        assert returncode == 0, f"expected success but mf6 failed:\n{out}"
    else:
        assert returncode != 0, f"expected failure but mf6 succeeded:\n{out}"
    for substring in required:
        assert substring in out, (
            f"expected message not found:\n  {substring!r}\nin output:\n{out}"
        )

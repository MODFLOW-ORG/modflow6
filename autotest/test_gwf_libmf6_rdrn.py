"""
Implements a "reverse drain" through the MODFLOW 6 API in a transient,
multi-stress-period simulation with pumping, and verifies the resulting heads and
flows for both the standard and the Newton-Raphson formulation.

A normal Drain (DRN) Package removes water from a cell when the head is above the
drain elevation. A reverse drain does the opposite: it injects water when the
head is below the elevation, driving the head up toward it.

A standard DRN Package cannot be reversed through the API, because the package
recomputes its own HCOF and RHS every formulate step. The reverse drain is
therefore implemented through the API Package, whose HCOF and RHS contributions
are supplied directly by the user and are not overwritten. Each outer (Picard or
Newton) iteration the head-dependent reverse-drain term

    q = COND * (ELEV - h)   when h < ELEV,   else 0

is written as HCOF = -COND and RHS = -COND * ELEV (or zero), under-relaxed
(damped) for stability.

Each (isolated) cell also has a well that extracts water. At steady state the
reverse-drain inflow balances the well outflow, and the head settles just below
the drain elevation, where the head-dependent inflow equals the extraction:

    h_steady = ELEV - |q_well| / COND

The head deficit below the elevation, |q_well| / COND, is inversely proportional
to the conductance: a head difference is required to drive the inflow, and a
larger COND drives the same inflow with a smaller difference. COND is therefore
chosen large here so the heads sit essentially on the drain elevations; a smaller
COND would hold them noticeably below, and COND -> infinity would pin them
exactly at the elevations (behaving like a constant-head boundary).

The well extraction rate is changed at the start of each stress period. Because
the reverse-drain inflow balances the well, the inflow tracks the pumping rate
from period to period, while the heads stay close to the (fixed) drain
elevations.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

NAME = "rdrn"
MODELNAME = NAME.upper()

# Five columns; columns 1 and 3 are inactive, so the three active cells
# (columns 0, 2, 4) are isolated and each reaches its own steady state.
idomain = [[[1, 0, 1, 0, 1]]]
active_cols = (0, 2, 4)
# reduced (one-based) node numbers of the active cells after the inactive
# columns are removed
active_nodes = np.array([1, 2, 3])

elevations = np.array([3.0, 5.0, 7.0])  # fixed reverse-drain elevations

# well extraction rate at each active cell, changed each stress period
well_rate_by_period = np.array([-1.0, -3.0, -2.0])
nper = well_rate_by_period.size

# reverse-drain conductance; the steady head sits |q_well| / cond below the
# elevation, so a large cond keeps the heads close to the drain elevations
cond = 1000.0
sy = 0.2
damping = 0.5  # under-relaxation of the reverse-drain term

# start at the period-1 steady-state heads so no storage is filled on step 1
strt = np.zeros((1, 1, 5))
strt[0, 0, list(active_cols)] = elevations - abs(well_rate_by_period[0]) / cond


def _build_model(ws, exe, newton):
    # exe is the libmf6 shared library so the test framework drives the run
    # through the API (api_func) rather than the mf6 executable
    sim = flopy.mf6.MFSimulation(
        sim_name=NAME, version="mf6", exe_name=str(exe), sim_ws=str(ws)
    )
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=[(20.0, 4, 1.0) for _ in range(nper)]
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=1e-9,
        outer_maximum=500,
        inner_dvclose=1e-11,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=NAME,
        newtonoptions="NEWTON" if newton else None,
        save_flows=True,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=5,
        delr=10.0,
        delc=10.0,
        top=10.0,
        botm=[0.0],
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0, print_flows=True)
    flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=0.0,
        sy=sy,
        transient={p: True for p in range(nper)},
    )
    # well extraction changes each stress period
    flopy.mf6.ModflowGwfwel(
        gwf,
        pname="wel",
        print_flows=True,
        stress_period_data={
            p: [[(0, 0, c), well_rate_by_period[p]] for c in active_cols]
            for p in range(nper)
        },
    )
    # the API Package supplies the reverse-drain HCOF and RHS at run time
    flopy.mf6.ModflowGwfapi(
        gwf, pname="rdrn", maxbound=active_nodes.size, print_flows=True
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{NAME}.hds",
        budget_filerecord=f"{NAME}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def _run_reverse_drain(libmf6, ws, results):
    """Drive the reverse drain through the API. Stash the final heads, per-step
    period, and inflow in ``results`` and return ``(success, stdout)`` for the
    test framework."""
    from modflowapi import ModflowApi

    output_file_path = os.path.join(str(ws), "mfsim.stdout")
    mf6 = ModflowApi(str(libmf6), working_directory=str(ws))
    mf6.initialize()

    nodelist = mf6.get_value_ptr(mf6.get_var_address("NODELIST", MODELNAME, "RDRN"))
    hcof = mf6.get_value_ptr(mf6.get_var_address("HCOF", MODELNAME, "RDRN"))
    rhs = mf6.get_value_ptr(mf6.get_var_address("RHS", MODELNAME, "RDRN"))
    nbound = mf6.get_value_ptr(mf6.get_var_address("NBOUND", MODELNAME, "RDRN"))
    head = mf6.get_value_ptr(mf6.get_var_address("X", MODELNAME))
    kper = mf6.get_value_ptr(mf6.get_var_address("KPER", "TDIS"))

    # assign the reverse drains to the active cells (constant for the run)
    nbound[0] = active_nodes.size
    nodelist[:] = active_nodes

    max_iter = int(mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))[0])

    period_per_step = []
    inflow_per_step = []
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        mf6.prepare_time_step(mf6.get_time_step())
        mf6.prepare_solve()
        kiter = 0
        while kiter < max_iter:
            # head-dependent reverse-drain term, evaluated at the current head
            h = head[active_nodes - 1]
            injecting = h < elevations
            new_hcof = np.where(injecting, -cond, 0.0)
            new_rhs = np.where(injecting, -cond * elevations, 0.0)
            # under-relax (damp) the term before each solve for stability
            hcof[:] = (1.0 - damping) * hcof + damping * new_hcof
            rhs[:] = (1.0 - damping) * rhs + damping * new_rhs
            if mf6.solve():
                break
            kiter += 1
        else:
            mf6.finalize()
            return False, open(output_file_path).readlines()
        h = head[active_nodes - 1]
        period_per_step.append(int(kper[0]))
        inflow_per_step.append(float(np.sum(cond * np.maximum(elevations - h, 0.0))))
        mf6.finalize_solve()
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

    results["heads"] = head[active_nodes - 1].copy()
    results["periods"] = np.array(period_per_step)
    results["inflow"] = np.array(inflow_per_step)
    mf6.finalize()
    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("newton", [False, True], ids=["standard", "newton"])
def test_reverse_drain(function_tmpdir, targets, newton):
    # the API run stashes its in-loop heads/periods/inflow here for check_output
    results = {}

    def build_models(test):
        return _build_model(test.workspace, test.targets["libmf6"], newton)

    def api_func(exe, ws):
        return _run_reverse_drain(exe, ws, results)

    def check_output(test):
        heads = results["heads"]
        periods = results["periods"]
        inflow = results["inflow"]

        # heads stay just below the (fixed) drain elevations
        expected_heads = elevations - abs(well_rate_by_period[-1]) / cond
        assert np.allclose(heads, expected_heads, atol=1e-3), (
            f"reverse-drain heads {heads} did not reach the expected steady-state "
            f"heads {expected_heads}"
        )

        # the total reverse-drain inflow tracks the well extraction in each period
        for p in range(nper):
            last_step = np.where(periods == p + 1)[0][-1]
            expected_inflow = abs(well_rate_by_period[p]) * active_nodes.size
            assert abs(inflow[last_step] - expected_inflow) < 1e-2, (
                f"period {p + 1}: reverse-drain inflow {inflow[last_step]} did not "
                f"balance the total well extraction {expected_inflow}"
            )

        # the inflow actually varied as the well rate changed
        assert inflow.max() - inflow.min() > 1.0, (
            "reverse-drain inflow did not vary with the changing well rate"
        )

    test = TestFramework(
        name=NAME,
        workspace=function_tmpdir,
        targets=targets,
        build=build_models,
        api_func=api_func,
        check=check_output,
        compare=None,
    )
    test.run()

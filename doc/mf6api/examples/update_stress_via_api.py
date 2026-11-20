"""
Updating stress-package input through the MODFLOW 6 API.

This example shows how to change a boundary input value at run time through the
API, and---more importantly---*when* the change must be applied. The same
pattern works whether the value is given as a constant in the package file or by
a time series; the difference is only in when MODFLOW 6 itself writes the value.

The model has a single General-Head Boundary (GHB) package with two boundaries:

  - boundary 0: conductance is a constant in the package file
  - boundary 1: conductance is supplied by a time series

MODFLOW 6 writes both values inside ``prepare_time_step``:

  - the package-file value is set once per stress period, during Read and Prepare
  - the time-series value is recomputed every time step, during Advance

Therefore an API override must be applied AFTER ``prepare_time_step`` (for
example, just before the solve). A value set before ``prepare_time_step`` is
overwritten---by Read and Prepare at the start of a stress period, or by the time
series on every step.

Run it with the paths to your MODFLOW 6 executable and shared library, e.g.::

    MF6_EXE=/path/to/mf6 MF6_LIB=/path/to/libmf6.so python update_stress_via_api.py
"""

import os
import tempfile

import flopy
import numpy as np
from modflowapi import ModflowApi

# Point these at your MODFLOW 6 executable and shared library.
MF6_EXE = os.environ.get("MF6_EXE", "mf6")
MF6_LIB = os.environ.get("MF6_LIB", "libmf6")

MODELNAME = "apiex"
GHB_PACKAGE = "ghb"
# the conductance we want to impose through the API, overriding both the
# package-file value (boundary 0) and the time-series value (boundary 1)
OVERRIDE_COND = 5.0


def build_model(ws):
    """Build a small transient model with a constant and a time-series GHB."""
    sim = flopy.mf6.MFSimulation(sim_name=MODELNAME, exe_name=MF6_EXE, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(10.0, 10, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-9,
        inner_dvclose=1e-10,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=MODELNAME, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=10.0, delc=10.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=8.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=1.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-4, sy=0.2, transient={0: True})
    # boundary 0: constant conductance (2.0) from the package file
    # boundary 1: conductance from the time series "condts" (ramps 1 -> 10)
    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        pname=GHB_PACKAGE,
        stress_period_data={
            0: [
                [(0, 0, 1), 6.0, 2.0],
                [(0, 0, 2), 6.0, "condts"],
            ]
        },
    )
    ghb.ts.initialize(
        filename=f"{GHB_PACKAGE}.ts",
        timeseries=[(0.0, 1.0), (10.0, 10.0)],
        time_series_namerecord=["condts"],
        interpolation_methodrecord=["linear"],
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 4), 1.0]]})
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{MODELNAME}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    sim.write_simulation(silent=True)


def run_with_override(ws):
    """Run the model, overriding the GHB conductance every time step."""
    mf6 = ModflowApi(MF6_LIB, working_directory=ws)
    mf6.initialize()

    # Pointer into the package's working conductance array (one entry per
    # boundary). Reading or writing through this pointer accesses live memory.
    cond = mf6.get_value_ptr(
        mf6.get_var_address("COND", MODELNAME.upper(), GHB_PACKAGE)
    )
    max_iter = int(mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))[0])

    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()
    while current_time < end_time:
        dt = mf6.get_time_step()

        # prepare_time_step performs Read and Prepare (start of a stress period)
        # and Advance (every step). On the first step of the stress period COND
        # holds the package-file value (boundary 0) and the time-series value
        # (boundary 1). On later steps the package-file value is left as it was
        # (so a previous override persists), while the time-series value is
        # recomputed and overwrites any previous override.
        mf6.prepare_time_step(dt)
        mf6.prepare_solve()

        # Inspect what MODFLOW 6 set, then override it. This MUST happen after
        # prepare_time_step; otherwise the time series (and Read and Prepare)
        # would overwrite our value. The time-series boundary must be overridden
        # every step; the package-file boundary only needs to be overridden once
        # per stress period, but overriding every step (as done here) is safe
        # for both.
        print(
            f"t={current_time:5.1f}  "
            f"file-based COND[0]={cond[0]:6.3f}  "
            f"time-series COND[1]={cond[1]:6.3f}  ->  "
            f"override both to {OVERRIDE_COND}"
        )
        cond[:] = OVERRIDE_COND

        kiter = 0
        while kiter < max_iter:
            if mf6.solve():
                break
            kiter += 1
        else:
            raise RuntimeError("did not converge")

        # the override is in effect for this solve, for BOTH boundaries
        assert np.allclose(cond, OVERRIDE_COND)

        mf6.finalize_solve()
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

    head = mf6.get_value(mf6.get_var_address("X", MODELNAME.upper())).copy()
    mf6.finalize()
    return head


def main():
    with tempfile.TemporaryDirectory() as ws:
        build_model(ws)
        head = run_with_override(ws)
        print("final head:", np.array2string(head, precision=3))


if __name__ == "__main__":
    main()

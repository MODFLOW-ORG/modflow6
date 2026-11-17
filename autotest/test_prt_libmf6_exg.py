"""
Run a GWF and PRT model in the same simulation through the API,
solving each solution in turn, and check that:

1. the particle tracks match a run of the same simulation with
   the mf6 executable
2. the committed particle store (PL* variables in each PRP's
   memory path), read after each time step, matches the track
   output up to the current time
3. pointers to the particle store remain valid across time steps
   in which no particles are released, i.e. while the store's
   size doesn't change

The flow system is the one on the FloPy readme. Particles are
released from the top left cell at two times, the second in
the middle of the simulation, so the store grows once.
"""

from pathlib import Path

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg
from prt_test_utils import FlopyReadmeCase, get_model_name

simname = "prtlibexg"
cases = [simname]
nstp = 5
releasetimes = [0.0, 0.5]  # second falls in the third time step
snapshot_file = "pl_snapshots.csv"

# particle store variables to check against the track output
pl_vars = {
    "PLIRPT": "irpt",
    "PLTRELEASE": "trelease",
    "PLTTRACK": "t",
    "PLX": "x",
    "PLY": "y",
    "PLZ": "z",
    "PLISTATUS": "istatus",
    "PLILAY": "ilay",
    "PLICU": "icell",
    "PLIZONE": "izone",
}


def build_sim(name, ws, exe):
    sim = FlopyReadmeCase.get_gwf_sim(name, ws, exe)
    sim.tdis.perioddata = [(FlopyReadmeCase.perlen, nstp, FlopyReadmeCase.tsmult)]

    prt_name = get_model_name(name, "prt")
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)
    flopy.mf6.ModflowGwfdis(
        prt,
        nlay=FlopyReadmeCase.nlay,
        nrow=FlopyReadmeCase.nrow,
        ncol=FlopyReadmeCase.ncol,
    )
    flopy.mf6.ModflowPrtmip(prt, porosity=FlopyReadmeCase.porosity)
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        nreleasepts=len(FlopyReadmeCase.releasepts_prt),
        packagedata=FlopyReadmeCase.releasepts_prt,
        nreleasetimes=len(releasetimes),
        releasetimes=[[t] for t in releasetimes],
        extend_tracking=True,
    )
    flopy.mf6.ModflowPrtoc(prt, trackcsv_filerecord=[f"{prt_name}.trk.csv"])
    flopy.mf6.ModflowGwfprt(
        sim,
        exgtype="GWF6-PRT6",
        exgmnamea=get_model_name(name, "gwf"),
        exgmnameb=prt_name,
    )
    ems = flopy.mf6.ModflowEms(sim, filename=f"{prt_name}.ems")
    sim.register_solution_package(ems, [prt.name])
    return sim


def build_models(idx, test):
    name = test.name
    api_sim = build_sim(name, test.workspace, test.targets["libmf6"])
    mf6_sim = build_sim(name, test.workspace / "mf6", test.targets["mf6"])
    return api_sim, mf6_sim


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx]
    model_ws = Path(model_ws or ".")
    prt_name = get_model_name(name, "prt").upper()
    output_file_path = model_ws / "mfsim.stdout"

    mf6 = ModflowApi(exe, working_directory=model_ws)
    mf6.initialize()

    def addr(var):
        return mf6.get_var_address(var, prt_name, "PRP1")

    nsol = mf6.get_subcomponent_count()
    assert nsol == 2, f"expected a GWF and a PRT solution, got {nsol}"
    max_iters = []
    for sol_id in range(1, nsol + 1):
        try:
            max_iters.append(
                mf6.get_value(mf6.get_var_address("MXITER", f"SLN_{sol_id}"))[0]
            )
        except Exception:
            # the explicit solution (PRT) has no MXITER, one solve is enough
            max_iters.append(1)

    snapshots = []
    ptrs = {}
    nptrchecks = 0
    kstp = 0
    try:
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.prepare_time_step(mf6.get_time_step())
            for sol_id, max_iter in enumerate(max_iters, start=1):
                mf6.prepare_solve(sol_id)
                for _ in range(max_iter):
                    if has_converged := mf6.solve(sol_id):
                        break
                mf6.finalize_solve(sol_id)
                if not has_converged:
                    return False, open(output_file_path).readlines()
            mf6.finalize_time_step()
            kstp += 1

            # pointers taken after the last step must still point to
            # the store if no particles were released this step
            new_ptrs = {var: mf6.get_value_ptr(addr(var)) for var in pl_vars}
            for var, ptr in ptrs.items():
                if ptr.size != new_ptrs[var].size:
                    continue
                assert np.shares_memory(ptr, new_ptrs[var]), (
                    f"step {kstp}: {var} was reallocated though its size is unchanged"
                )
                nptrchecks += 1
            ptrs = new_ptrs

            snapshot = pd.DataFrame(
                {col: mf6.get_value(addr(var)) for var, col in pl_vars.items()}
            )
            snapshot["kstp"] = kstp
            snapshot["time"] = mf6.get_current_time()
            snapshots.append(snapshot)
    finally:
        mf6.finalize()

    assert kstp == nstp
    # releases on steps 1 and 3 only, so the store's size holds on 2, 4 and 5
    assert nptrchecks == 3 * len(pl_vars)
    pd.concat(snapshots).to_csv(model_ws / snapshot_file, index=False)
    return True, open(output_file_path).readlines()


def check_output(idx, test):
    name = test.name
    prt_name = get_model_name(name, "prt")
    trk_name = f"{prt_name}.trk.csv"
    api_ws = test.workspace
    mf6_ws = test.workspace / "mf6"

    # tracks from the API run match tracks from the mf6 run
    api_trk = pd.read_csv(api_ws / trk_name)
    mf6_trk = pd.read_csv(mf6_ws / trk_name)
    numeric = api_trk.select_dtypes("number").columns
    assert api_trk.shape == mf6_trk.shape
    assert np.allclose(api_trk[numeric], mf6_trk[numeric], equal_nan=True)

    # every release is in the tracks
    nrpt = len(FlopyReadmeCase.releasepts_prt)
    pids = ["iprp", "irpt", "trelease"]
    assert len(api_trk.groupby(pids)) == nrpt * len(releasetimes)

    # the particle store matches the last track record at or before
    # the end of each time step. on the last step, extended tracking
    # carries particles past the end, so use their last record.
    snapshots = pd.read_csv(api_ws / snapshot_file)
    assert sorted(snapshots.kstp.unique()) == list(range(1, nstp + 1))
    for kstp, snapshot in snapshots.groupby("kstp"):
        tcur = snapshot.time.iloc[0]
        released = [t for t in releasetimes if t <= tcur]
        assert len(snapshot) == nrpt * len(released), (
            f"step {kstp}: expected {nrpt * len(released)} particles, "
            f"got {len(snapshot)}"
        )
        upto = api_trk if kstp == nstp else api_trk[api_trk.t <= tcur * (1 + 1e-12)]
        expected = upto.groupby(["irpt", "trelease"], as_index=False).last()
        actual = snapshot.sort_values(["irpt", "trelease"]).reset_index(drop=True)
        expected = expected.sort_values(["irpt", "trelease"]).reset_index(drop=True)
        assert len(actual) == len(expected), f"step {kstp}"
        for col in pl_vars.values():
            assert np.allclose(actual[col], expected[col]), (
                f"step {kstp}: {col} differs from track output\n"
                f"store:\n{actual[col].to_numpy()}\n"
                f"tracks:\n{expected[col].to_numpy()}"
            )


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()

"""
Verify how static (PACKAGEDATA and OUTLETS) values written through the API
interact with time series in the advanced GWF packages.

The packages copy time series capable static values from the input into
their own arrays each time step. Only values linked to a time series may be
copied; any other value written through the API must persist.

Each case writes two entries of one variable, in the same input column:

    case   package  written (old name)  entry L       entry U
    maw    MAW      WELL_HEAD           well 2 STRT   well 1 STRT
    sfr    SFR      HK                  reach 2 RHK   reach 1 RHK
    lak    LAK      OUTWIDTH            outlet 1      outlet 2

Entry U is never linked to a time series. Entry L is one of:

    kind      entry L                time series file
    plain     numeric value          none
    plain_ts  numeric value          drives a PERIOD entry
    flat      linked, one value      drives entry L
    ramp      linked, every step     drives entry L

MAW wells are constant head. PACKAGEDATA and OUTLETS rows are listed out of
feature order, so the input row of each entry differs from its feature
number. The write happens once, during the first time step, at one of three
points:

    phase   written after
    init    initialize
    pts     prepare_time_step
    ps      prepare_solve

Time series advance at the end of prepare_time_step, and the packages copy
linked values into their own arrays in prepare_solve, so a linked entry
written after init or pts is overwritten within the first time step, and one
written after ps holds for the first time step only. An entry that is not
linked keeps the written value, and the solution matches a reference run
whose input sets that value.
"""

import flopy
import numpy as np
import pytest
from modflow_devtools.markers import requires_pkg

cases = ["maw", "sfr", "lak"]
kinds = ["plain", "plain_ts", "flat", "ramp"]
phases = ["init", "pts", "ps"]

name = "api_st_writes"
nper, perlen, nstp = 2, 4.0, 2
top, botm = 10.0, 0.0
tsname = "linked"

# -- case: (package, written old name, package state to compare,
#           entry L index, written value, flat series value,
#           ramp series values at each step start,
#           entry U index, input value, written value,
#           PERIOD entry value driven by the plain_ts series)
writes = {
    "maw": (
        "MAW-1",
        "WELL_HEAD",
        "HEAD",
        1,
        7.0,
        9.0,
        [9.0, 8.5, 8.0, 7.5],
        0,
        9.0,
        8.0,
        -50.0,
    ),
    "sfr": (
        "SFR-1",
        "HK",
        "DSFLOW",
        1,
        0.0,
        0.1,
        [0.1, 0.2, 0.3, 0.4],
        0,
        0.1,
        0.05,
        10.0,
    ),
    "lak": (
        "LAK-1",
        "OUTWIDTH",
        "SIMOUTRATE",
        0,
        10.0,
        2.0,
        [2.0, 3.0, 4.0, 5.0],
        1,
        1.0,
        5.0,
        0.5,
    ),
}


def linked(kind):
    return kind in ("flat", "ramp")


def input_values(case, kind):
    """Value of entry L at each time step without a write."""
    _, _, _, _, _, flat, ramp, _, _, _, _ = writes[case]
    return ramp if kind == "ramp" else [flat] * (nper * nstp)


def get_sim(ws, case, kind, reference):
    """Build the simulation. For the reference run, the input sets each entry
    that is not linked to the value the API run writes."""
    _, _, _, _, lvalue, flat, _, _, uinput, uvalue, other = writes[case]
    lin = tsname if linked(kind) else (lvalue if reference else flat)
    uin = uvalue if reference else uinput
    pin = tsname if kind == "plain_ts" else other

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=[(perlen, nstp, 1.0)] * nper)
    flopy.mf6.ModflowIms(
        sim,
        complexity="complex",
        outer_dvclose=1e-9,
        inner_dvclose=1e-9,
        outer_maximum=200,
        inner_maximum=200,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    ncol = 5 if case == "maw" else 3
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=ncol, delr=100.0, delc=100.0, top=top, botm=botm
    )
    k = 10.0 if case == "maw" else 1.0
    strt = 9.0 if case == "maw" else 8.0
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)

    if case == "maw":
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: [[(0, 0, 0), strt], [(0, 0, ncol - 1), strt]]},
        )
        pd = [[1, 0.1, botm, lin, "THIEM", 1], [0, 0.1, botm, uin, "THIEM", 1]]
        cd = [[n, 0, (0, 0, 1 + 2 * n), top, botm, -999.0, -999.0] for n in range(2)]
        per = [[0, "STATUS", "CONSTANT"], [1, "STATUS", "CONSTANT"]]
        pkg = flopy.mf6.ModflowGwfmaw(
            gwf,
            nmawwells=2,
            packagedata=pd,
            connectiondata=cd,
            perioddata={0: per, 1: [[0, "RATE", pin]]},
            pname="MAW-1",
        )

    elif case == "sfr":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 0), strt]]})
        rhk = [uin, lin, 0.1]
        ncon = [1, 2, 1]
        pd = [
            [i, (0, 0, i), 100.0, 5.0, 0.001, 9.0, 1.0, rhk[i], 0.03, ncon[i], 1.0, 0]
            for i in [1, 2, 0]
        ]
        cd = [[0, -1], [1, 0, -2], [2, 1]]
        per = [[0, "INFLOW", pin]]
        pkg = flopy.mf6.ModflowGwfsfr(
            gwf,
            nreaches=3,
            packagedata=pd,
            connectiondata=cd,
            perioddata={0: per, 1: per},
            pname="SFR-1",
        )

    elif case == "lak":
        flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 2), strt]]})
        pd = [[0, 11.0, 1]]
        cd = [[0, 0, (0, 0, 0), "VERTICAL", 0.001, 0.0, 0.0, 0.0, 0.0]]
        outlets = [
            [1, 0, -1, "WEIR", 10.8, uin, 0.03, 0.001],
            [0, 0, -1, "WEIR", 10.5, lin, 0.03, 0.001],
        ]
        per = [[0, "RAINFALL", pin]]
        pkg = flopy.mf6.ModflowGwflak(
            gwf,
            nlakes=1,
            noutlets=2,
            packagedata=pd,
            connectiondata=cd,
            outlets=outlets,
            perioddata={0: per, 1: per},
            pname="LAK-1",
        )

    if kind != "plain":
        # -- stepwise series with a breakpoint at the start of each time step
        delt = perlen / nstp
        if kind == "plain_ts":
            values = [other] * (nper * nstp)
        else:
            values = input_values(case, kind)
        timeseries = [(i * delt, v) for i, v in enumerate(values)]
        timeseries.append((nper * perlen, values[-1]))
        pkg.ts.initialize(
            filename=f"{name}.{case}.ts",
            timeseries=timeseries,
            time_series_namerecord=tsname,
            interpolation_methodrecord="stepwise",
        )
    return sim


def run_sim(lib, ws, case, phase=None):
    """Run through the API, writing entries L and U once at the given phase of
    the first time step (no write if phase is None). Returns per-step model
    heads, package state, and entries L and U."""
    from xmipy import XmiWrapper

    pname, var, state, lidx, lvalue, _, _, uidx, _, uvalue, _ = writes[case]
    prefix = f"{name.upper()}/{pname}"
    mf6 = XmiWrapper(lib, working_directory=str(ws))
    mf6.initialize()
    kstp = 0

    def write(at):
        if at == phase and kstp == 0:
            ptr = mf6.get_value_ptr(mf6.get_var_address(var, name.upper(), pname))
            ptr[lidx] = lvalue
            ptr[uidx] = uvalue

    try:
        write("init")
        heads, states, entries = [], [], []
        while mf6.get_current_time() < mf6.get_end_time():
            mf6.prepare_time_step(mf6.get_time_step())
            write("pts")
            nsol = mf6.get_subcomponent_count()
            for sol_id in range(1, nsol + 1):
                mf6.prepare_solve(sol_id)
                if sol_id == nsol:
                    write("ps")
                for _ in range(mf6.get_value("SLN_1/MXITER")[0]):
                    if mf6.solve(sol_id):
                        break
                else:
                    raise RuntimeError(f"no convergence at step {kstp + 1}")
                mf6.finalize_solve(sol_id)
            mf6.finalize_time_step()
            kstp += 1
            heads.append(mf6.get_value(f"{name.upper()}/X"))
            states.append(mf6.get_value(f"{prefix}/{state}"))
            v = mf6.get_value(f"{prefix}/{var}")
            entries.append((v[lidx], v[uidx]))
    finally:
        mf6.finalize()
    entries = np.array(entries)
    return np.array(heads), np.array(states), entries[:, 0], entries[:, 1]


def check(lib, ws, case, kind, phase):
    _, var, state, lidx, lvalue, _, _, uidx, _, uvalue, _ = writes[case]
    ref_ws, api_ws = ws / "reference", ws / "api"
    get_sim(ref_ws, case, kind, reference=True).write_simulation(silent=True)
    get_sim(api_ws, case, kind, reference=False).write_simulation(silent=True)
    ref_heads, ref_states, ref_l, _ = run_sim(lib, ref_ws, case)
    heads, states, written_l, written_u = run_sim(lib, api_ws, case, phase)

    # -- a linked entry L follows its time series, except that a write after
    #    prepare_solve holds for the first time step. Otherwise entry L keeps
    #    the written value
    if linked(kind):
        values = np.array(input_values(case, kind))
        assert np.allclose(ref_l, values), (
            f"{var}[{lidx}] does not follow the time series: {ref_l} (series {values})"
        )
        expected = values.copy()
        if phase == "ps":
            expected[0] = lvalue
    else:
        expected = np.full(nper * nstp, lvalue)
    assert np.allclose(written_l, expected), (
        f"{var}[{lidx}] written after {phase} ({kind}): {written_l} after "
        f"each step (expected {expected})"
    )

    # -- entry U is never linked and keeps the written value
    assert np.allclose(written_u, uvalue), (
        f"unlinked {var}[{uidx}] written after {phase} ({kind}) was "
        f"overwritten: {written_u} after each step (wrote {uvalue})"
    )

    # -- the written values reach the solution. The reference run can't hold a
    #    linked entry at the written value, so skip the first time step when
    #    that entry was written after prepare_solve
    first = 1 if linked(kind) and phase == "ps" else 0
    assert np.allclose(states[first:], ref_states[first:]), (
        f"{state} differs from the reference run: {states} vs {ref_states}"
    )
    assert np.allclose(heads[first:], ref_heads[first:]), (
        f"heads differ from the reference run: {heads} vs {ref_heads}"
    )


@requires_pkg("xmipy")
@pytest.mark.parametrize("phase", phases)
@pytest.mark.parametrize("kind", kinds)
@pytest.mark.parametrize("case", cases)
def test_api_static_writes(case, kind, phase, function_tmpdir, targets):
    check(targets["libmf6"], function_tmpdir, case, kind, phase)

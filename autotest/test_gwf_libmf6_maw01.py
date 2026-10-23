"""
Verify that deprecated MAW variable names remain resolvable through the API.

IDM integration moved several MAW well variables from the package memory path
into the input context, renaming one in the process. MAW registers aliases at
the package path so external programs (modflowapi, coupled models) can still
address these variables by their old names.

    old name   ->  current name (input context)
    RATE       ->  RATE        (path-only move)
    WELL_HEAD  ->  WELL_HEAD   (path-only move)
    MAUXVAR    ->  AUX         (renamed)

For each: the old name resolves, its value matches the model input, and it
shares memory with the current input-context variable. Each old name used is
reported once as deprecated in the simulation listing file.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_maw01"]

# -- MAW package name (as seen in the memory manager / API)
maw_pname = "MAW-1"

# -- deprecated names read through the API below
used_old_names = ["RATE", "WELL_HEAD", "MAUXVAR"]

# -- grid
nlay, nrow, ncol = 1, 1, 5
delr = delc = 100.0
top = 10.0
botm = 0.0
strt = 9.0
k11 = 10.0

# -- well input values (distinct so each alias can be told apart)
nwells = 2
radius = 0.1
rates = [-50.0, 0.0]  # -> RATE; well 2 is constant head
well_head = 8.5  # -> WELL_HEAD for well 2
aux_vals = [10.0, 11.0]  # -> AUX (old MAUXVAR)

nouter, ninner = 100, 50
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def get_model(ws, name):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="all",
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
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
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=k11)
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(0, 0, 0), strt], [(0, 0, ncol - 1), strt]]}
    )

    # -- MAW: one connection per well, in columns 2 and 4
    packagedata = [
        [n, radius, botm, strt, "THIEM", 1, aux_vals[n]] for n in range(nwells)
    ]
    connectiondata = [
        [n, 0, (0, 0, 1 + 2 * n), top, botm, -999.0, -999.0] for n in range(nwells)
    ]
    perioddata = {
        0: [
            [0, "RATE", rates[0]],
            [1, "STATUS", "CONSTANT"],
            [1, "WELL_HEAD", well_head],
        ]
    }
    flopy.mf6.ModflowGwfmaw(
        gwf,
        auxiliary=["AUX1"],
        nmawwells=nwells,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname=maw_pname,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sim = get_model(test.workspace, name)
    # comparison sim in a "libmf6" subdir: the framework runs this one via the API
    sim_api = get_model(os.path.join(test.workspace, "libmf6"), name)
    return sim, sim_api


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx].upper()
    if model_ws is None:
        model_ws = "."
    output_file_path = os.path.join(model_ws, "mfsim.stdout")

    try:
        mf6 = ModflowApi(exe, working_directory=model_ws)
    except Exception as e:
        print("Failed to load " + str(exe))
        print("with message: " + str(e))
        return False, open(output_file_path).readlines()

    try:
        mf6.initialize()
    except Exception:
        return False, open(output_file_path).readlines()

    try:
        # advance one step so period (RATE/WELL_HEAD) data is populated
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()
        mxit_tag = mf6.get_var_address("MXITER", "SLN_1")
        max_iter = mf6.get_value(mxit_tag)
        while current_time < end_time:
            mf6.prepare_time_step(mf6.get_time_step())
            mf6.prepare_solve()
            kiter = 0
            while kiter < max_iter:
                has_converged = mf6.solve()
                kiter += 1
                if has_converged:
                    break
            if not has_converged:
                print("model did not converge")
                return False, open(output_file_path).readlines()
            mf6.finalize_solve()
            mf6.finalize_time_step()
            current_time = mf6.get_current_time()

        input_prefix = f"__INPUT__/{name}/{maw_pname}"

        # -- RATE (moved): value and shared memory with the input context
        rate_tag = mf6.get_var_address("RATE", name, maw_pname)
        rate = mf6.get_value(rate_tag)
        assert np.allclose(rate, rates), f"old name 'RATE' value {rate} != {rates}"
        assert np.allclose(rate, mf6.get_value(f"{input_prefix}/RATE")), (
            "alias 'RATE' differs from input-context 'RATE'"
        )

        # -- prove shared storage: write through the old name, read back the
        #    input-context variable
        rate_mod = rate.copy()
        rate_mod[0] -= 1.0
        mf6.set_value(rate_tag, rate_mod)
        rate_after = mf6.get_value(f"{input_prefix}/RATE")
        assert np.isclose(rate_after[0], rates[0] - 1.0), (
            "writing through deprecated name 'RATE' was not reflected in "
            f"input-context 'RATE' (got {rate_after[0]})"
        )
        rate_mod[0] = rates[0]
        mf6.set_value(rate_tag, rate_mod)

        # -- WELL_HEAD (moved): constant-head value for well 2
        whead = mf6.get_value(mf6.get_var_address("WELL_HEAD", name, maw_pname))
        assert np.isclose(whead[1], well_head), (
            f"old name 'WELL_HEAD' well-2 value {whead[1]} != {well_head}"
        )
        assert np.allclose(whead, mf6.get_value(f"{input_prefix}/WELL_HEAD")), (
            "alias 'WELL_HEAD' differs from input-context 'WELL_HEAD'"
        )

        # -- well auxiliary values: MAUXVAR -> AUX (renamed)
        mauxvar = mf6.get_value(mf6.get_var_address("MAUXVAR", name, maw_pname))
        assert np.allclose(np.ravel(mauxvar), aux_vals), (
            f"old name 'MAUXVAR' value {np.ravel(mauxvar)} != {aux_vals}"
        )
        aux = mf6.get_value(f"{input_prefix}/AUX")
        assert np.allclose(mauxvar, aux), (
            f"alias 'MAUXVAR' and 'AUX' differ: {mauxvar} vs {aux}"
        )

    except Exception as e:
        print("Deprecated-name resolution check failed:")
        print(f"{type(e).__name__}: {e}")
        return False, open(output_file_path).readlines()

    try:
        mf6.finalize()
    except Exception:
        return False, open(output_file_path).readlines()

    # -- each old name used above is reported once
    expected = {f"{old}@{name}/{maw_pname}": 1 for old in used_old_names}
    reported = deprecation_warnings(os.path.join(model_ws, "mfsim.lst"))
    if reported != expected:
        print(f"deprecation warnings {reported} != expected {expected}")
        return False, open(output_file_path).readlines()

    return True, open(output_file_path).readlines()


def deprecation_warnings(lst_path):
    """Count deprecated-name warnings in a listing file, keyed by name@path."""
    # warnings wrap across lines, so match on whitespace-collapsed text
    with open(lst_path) as f:
        text = " ".join(f.read().split())
    counts = {}
    pattern = r"API VARIABLE '(\S+)' IN '(\S+)' WAS DEPRECATED"
    for var, path in re.findall(pattern, text):
        key = f"{var}@{path}"
        counts[key] = counts.get(key, 0) + 1
    return counts


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()

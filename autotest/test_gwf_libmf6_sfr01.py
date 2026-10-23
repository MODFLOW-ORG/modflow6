"""
Verify that deprecated SFR variable names remain resolvable through the API.

IDM integration moved several SFR reach variables from the package memory path
into the input context, renaming some in the process. SFR registers aliases at
the package path so external programs (modflowapi, coupled models) can still
address these variables by their old names.

A representative subset is checked:

    old name  ->  current name (input context)
    HK        ->  BEDK            (renamed)
    ROUGH     ->  MANNING         (renamed)
    USTRF     ->  UPSTREAM_FRAC   (renamed)
    INFLOW    ->  INFLOW          (path-only move)
    RAIN      ->  RAINFALL        (renamed)
    RAUXVAR   ->  AUX             (renamed)

For each: the old name resolves, its value matches the model input, and it
shares memory with the current input-context variable. Each old name used is
reported once as deprecated in the simulation listing file; aliases that are
not used (EVAP, RUNOFF, SSTAGE) are not reported. Without diversions DIVFLOW
still resolves, as an empty package array that is not reported.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_sfr01"]

# -- SFR package name (as seen in the memory manager / API)
sfr_pname = "SFR-1"

# -- deprecated names read through the API below
used_old_names = ["HK", "ROUGH", "USTRF", "INFLOW", "RAIN", "RAUXVAR"]

# -- grid
nlay, nrow, ncol = 1, 1, 5
delr = delc = 100.0
top = 10.0
botm = 0.0
strt = 9.0
k11 = 1.0

# -- reach input values (distinct so each alias can be told apart)
nreaches = ncol
rlen = delr
rwid = 5.0
rgrd = 0.001
rtp = 5.0
rbth = 1.0
rhk = 0.1  # -> BEDK (old HK)
man = 0.035  # -> MANNING (old ROUGH)
ustrf = 1.0  # -> UPSTREAM_FRAC (old USTRF)
aux_vals = [10.0 + irch for irch in range(nreaches)]  # -> AUX (old RAUXVAR)

# -- period input values
inflow_rate = 100.0  # -> INFLOW (old INFLOW)
rain_rate = 0.5  # -> RAINFALL (old RAIN)

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

    # -- SFR: one reach per column, linear chain
    packagedata = []
    for irch in range(nreaches):
        nconn = 1 if irch in (0, nreaches - 1) else 2
        packagedata.append(
            [
                irch,
                (0, 0, irch),
                rlen,
                rwid,
                rgrd,
                rtp,
                rbth,
                rhk,
                man,
                nconn,
                ustrf,
                0,
                aux_vals[irch],
            ]
        )
    connectiondata = []
    for irch in range(nreaches):
        rc = [irch]
        if irch > 0:
            rc.append(irch - 1)
        if irch < nreaches - 1:
            rc.append(-(irch + 1))
        connectiondata.append(rc)

    perioddata = {
        0: [
            [0, "INFLOW", inflow_rate],
            *[[irch, "RAINFALL", rain_rate] for irch in range(nreaches)],
        ]
    }

    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        auxiliary=["AUX1"],
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname=sfr_pname,
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

    # advance one step so period (RAINFALL/INFLOW) data is populated
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()

    try:
        # (old name, current name, expected scalar value)
        checks = [
            ("HK", "BEDK", rhk),
            ("ROUGH", "MANNING", man),
            ("USTRF", "UPSTREAM_FRAC", ustrf),
        ]

        while current_time < end_time:
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0
            mf6.prepare_solve()
            mxit_tag = mf6.get_var_address("MXITER", "SLN_1")
            max_iter = mf6.get_value(mxit_tag)
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

        # -- old names resolve, carry the input value, and share memory with
        #    the current input-context variable
        input_prefix = f"__INPUT__/{name}/{sfr_pname}"
        for old_name, new_name, expected in checks:
            old_tag = mf6.get_var_address(old_name, name, sfr_pname)
            old_val = mf6.get_value(old_tag)
            assert old_val is not None and len(old_val) == nreaches, (
                f"old name '{old_name}' did not resolve to a reach array"
            )
            assert np.allclose(old_val, expected), (
                f"old name '{old_name}' value {old_val} != expected {expected}"
            )

            # -- alias shares memory with the input-context variable
            new_val = mf6.get_value(f"{input_prefix}/{new_name}")
            assert np.allclose(old_val, new_val), (
                f"alias '{old_name}' and '{new_name}' differ: {old_val} vs {new_val}"
            )

        # -- prove shared storage: write through the old name, read back the
        #    input-context variable
        hk_tag = mf6.get_var_address("HK", name, sfr_pname)
        hk_mod = mf6.get_value(hk_tag).copy()
        hk_mod[0] += 1.0
        mf6.set_value(hk_tag, hk_mod)
        bedk_after = mf6.get_value(f"{input_prefix}/BEDK")
        assert np.isclose(bedk_after[0], rhk + 1.0), (
            "writing through deprecated name 'HK' was not reflected in "
            f"input-context 'BEDK' (got {bedk_after[0]}, expected {rhk + 1.0})"
        )
        # -- restore
        hk_mod[0] = rhk
        mf6.set_value(hk_tag, hk_mod)

        # -- period variables, populated after the step: INFLOW (moved),
        #    RAIN -> RAINFALL (renamed)
        inflow_tag = mf6.get_var_address("INFLOW", name, sfr_pname)
        inflow_val = mf6.get_value(inflow_tag)
        assert np.isclose(inflow_val[0], inflow_rate), (
            f"old name 'INFLOW' reach-1 value {inflow_val[0]} != {inflow_rate}"
        )

        rain_tag = mf6.get_var_address("RAIN", name, sfr_pname)
        rain_val = mf6.get_value(rain_tag)
        assert np.allclose(rain_val, rain_rate), (
            f"old name 'RAIN' value {rain_val} != {rain_rate}"
        )
        rainfall_val = mf6.get_value(f"{input_prefix}/RAINFALL")
        assert np.allclose(rain_val, rainfall_val), (
            f"alias 'RAIN' and 'RAINFALL' differ: {rain_val} vs {rainfall_val}"
        )

        # -- reach auxiliary values: RAUXVAR -> AUX (renamed)
        rauxvar = mf6.get_value(mf6.get_var_address("RAUXVAR", name, sfr_pname))
        assert np.allclose(np.ravel(rauxvar), aux_vals), (
            f"old name 'RAUXVAR' value {np.ravel(rauxvar)} != {aux_vals}"
        )
        aux = mf6.get_value(f"{input_prefix}/AUX")
        assert np.allclose(rauxvar, aux), (
            f"alias 'RAUXVAR' and 'AUX' differ: {rauxvar} vs {aux}"
        )

        # -- no diversions: DIVFLOW is an empty package array, as before IDM
        divflow_tag = mf6.get_var_address("DIVFLOW", name, sfr_pname)
        assert mf6.get_var_nbytes(divflow_tag) == 0, (
            "old name 'DIVFLOW' is not an empty array without diversions"
        )

    except Exception as e:
        print("Deprecated-name resolution check failed:")
        print(f"{type(e).__name__}: {e}")
        return False, open(output_file_path).readlines()

    try:
        mf6.finalize()
    except Exception:
        return False, open(output_file_path).readlines()

    # -- each old name used above is reported once; unused aliases are not
    expected = {f"{old}@{name}/{sfr_pname}": 1 for old in used_old_names}
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

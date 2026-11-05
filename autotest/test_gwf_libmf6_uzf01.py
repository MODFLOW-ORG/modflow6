"""
Verify that deprecated UZF variable names remain resolvable through the API.

IDM integration moved the UZF time-series aware cell variables and auxiliary
array from the package memory path into the input context, renaming them. UZF
registers aliases at the package path so external programs (modflowapi,
coupled models) can still address these variables by their old names.

    old name      ->  current name (input context)
    SINF_PVAR     ->  FINF
    PET_PVAR      ->  PET
    EXDP_PVAR     ->  EXTDP
    EXTWC_PVAR    ->  EXTWC
    HA_PVAR       ->  HA
    HROOT_PVAR    ->  HROOT
    ROOTACT_PVAR  ->  ROOTACT
    UAUXVAR       ->  AUXVAR

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

cases = ["libgwf_uzf01"]

# -- UZF package name (as seen in the memory manager / API)
uzf_pname = "UZF-1"

# -- (old name, current name, period data value)
renames = [
    ("SINF_PVAR", "FINF", 0.01),
    ("PET_PVAR", "PET", 0.002),
    ("EXDP_PVAR", "EXTDP", 1.5),
    ("EXTWC_PVAR", "EXTWC", 0.06),
    ("HA_PVAR", "HA", -10.0),
    ("HROOT_PVAR", "HROOT", -20.0),
    ("ROOTACT_PVAR", "ROOTACT", 0.5),
]
aux_val = 7.0  # -> AUXVAR (old UAUXVAR)

# -- grid
nlay, nrow, ncol = 1, 1, 3
delr = delc = 100.0
top = 10.0
botm = 0.0
strt = 5.0
k11 = 1.0

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
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=k11)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient={0: True})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(0, 0, 0), strt]]})

    # -- UZF: one cell per column, each on the land surface; auxiliary values
    #    are period data in UZF
    packagedata = [
        [n, (0, 0, n), 1, -1, 0.1, k11, 0.05, 0.2, 0.05, 4.0] for n in range(ncol)
    ]
    perioddata = {
        0: [[n] + [v for _, _, v in renames] + [aux_val] for n in range(ncol)]
    }
    flopy.mf6.ModflowGwfuzf(
        gwf,
        auxiliary=["AUX1"],
        simulate_et=True,
        ntrailwaves=7,
        nwavesets=40,
        nuzfcells=ncol,
        packagedata=packagedata,
        perioddata=perioddata,
        pname=uzf_pname,
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
        # advance one step so period data is populated
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

        input_prefix = f"__INPUT__/{name}/{uzf_pname}"

        # -- old names carry the period input and share memory with the
        #    input-context variable
        for old_name, new_name, value in renames:
            old_val = mf6.get_value(mf6.get_var_address(old_name, name, uzf_pname))
            assert np.allclose(old_val, value), (
                f"old name '{old_name}' value {old_val} != {value}"
            )
            new_val = mf6.get_value(f"{input_prefix}/{new_name}")
            assert np.array_equal(old_val, new_val), (
                f"alias '{old_name}' and '{new_name}' differ: {old_val} vs {new_val}"
            )

        # -- cell auxiliary values: UAUXVAR -> AUXVAR (renamed)
        uauxvar = mf6.get_value(mf6.get_var_address("UAUXVAR", name, uzf_pname))
        assert np.allclose(uauxvar, aux_val), (
            f"old name 'UAUXVAR' value {uauxvar} != {aux_val}"
        )
        assert np.array_equal(uauxvar, mf6.get_value(f"{input_prefix}/AUXVAR")), (
            "alias 'UAUXVAR' differs from input-context 'AUXVAR'"
        )

        # -- prove shared storage: write through the old name, read back the
        #    input-context variable
        tag = mf6.get_var_address("SINF_PVAR", name, uzf_pname)
        finf = mf6.get_value(tag)
        mf6.set_value(tag, finf + 1.0)
        after = mf6.get_value(f"{input_prefix}/FINF")
        assert np.allclose(after, finf + 1.0), (
            "writing through deprecated name 'SINF_PVAR' was not reflected in "
            "input-context 'FINF'"
        )
        mf6.set_value(tag, finf)

    except Exception as e:
        print("Deprecated-name resolution check failed:")
        print(f"{type(e).__name__}: {e}")
        return False, open(output_file_path).readlines()

    try:
        mf6.finalize()
    except Exception:
        return False, open(output_file_path).readlines()

    # -- each old name used above is reported once
    used = [old for old, _, _ in renames] + ["UAUXVAR"]
    expected = {f"{old}@{name}/{uzf_pname}": 1 for old in used}
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

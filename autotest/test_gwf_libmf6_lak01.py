"""
Verify that deprecated LAK variable names remain resolvable through the API.

IDM integration moved the LAK lake and outlet period variables from the package
memory path into the input context, renaming the outlet variables and the
auxiliary array, and renamed the lake connection cell array. The legacy solver
refactor renamed the fallback variables. LAK registers aliases at the package
path so external programs (modflowapi, coupled models) can still address these
variables by their old names. All aliases are checked:

    old name    ->  current name
    STAGE       ->  STAGE      (input context, path-only move)
    RAINFALL    ->  RAINFALL   (input context, path-only move)
    EVAPORATION ->  EVAPORATION (input context, path-only move)
    RUNOFF      ->  RUNOFF     (input context, path-only move)
    INFLOW      ->  INFLOW     (input context, path-only move)
    WITHDRAWAL  ->  WITHDRAWAL (input context, path-only move)
    OUTRATE     ->  RATE       (input context, renamed)
    OUTINVERT   ->  INVERT     (input context, renamed)
    OUTWIDTH    ->  WIDTH      (input context, renamed)
    OUTROUGH    ->  ROUGH      (input context, renamed)
    OUTSLOPE    ->  SLOPE      (input context, renamed)
    LAUXVAR     ->  AUX        (input context, renamed)
    CELLID      ->  GWFCONN    (package path, renamed)
    IFALLBACK   ->  ILEGACY    (package path, renamed)
    IFORCEFB    ->  IFORCELEG  (package path, renamed)

For each: the old name resolves, its value matches the model input, and it
shares memory with the current variable. Each old name used is reported once
as deprecated in the simulation listing file. A second case without outlets or
auxiliary variables checks that, as in 6.8.1, LAUXVAR resolves as an empty
array and the outlet names do not exist.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_lak01", "libgwf_lak01b"]

# -- case 1 has an outlet and an auxiliary variable; case 2 has neither
with_extras = [True, False]

# -- LAK package name (as seen in the memory manager / API)
lak_pname = "LAK-1"

# -- outlet names, only present when outlets are
outlet_names = ["OUTRATE", "OUTINVERT", "OUTWIDTH", "OUTROUGH", "OUTSLOPE"]

# -- (old name, current name, current name is in the input context)
renames = [
    ("STAGE", "STAGE", True),
    ("RAINFALL", "RAINFALL", True),
    ("EVAPORATION", "EVAPORATION", True),
    ("RUNOFF", "RUNOFF", True),
    ("INFLOW", "INFLOW", True),
    ("WITHDRAWAL", "WITHDRAWAL", True),
    ("OUTRATE", "RATE", True),
    ("OUTINVERT", "INVERT", True),
    ("OUTWIDTH", "WIDTH", True),
    ("OUTROUGH", "ROUGH", True),
    ("OUTSLOPE", "SLOPE", True),
    ("LAUXVAR", "AUX", True),
    ("CELLID", "GWFCONN", False),
    ("IFALLBACK", "ILEGACY", False),
    ("IFORCEFB", "IFORCELEG", False),
]

# -- grid
nlay, nrow, ncol = 1, 1, 5
delr = delc = 100.0
top = 10.0
botm = 0.0
strt = 9.0
k11 = 1.0

# -- lake and outlet input values (distinct so each alias can be told apart)
lake_strt = 9.0
lake_col = 2  # lake connects horizontally to this column
rainfall = 0.5
evaporation = 0.25
runoff = 3.0
inflow = 4.0
withdrawal = 1.5
aux_val = 7.0  # -> AUX (old LAUXVAR)
invert, width, rough, slope = 8.5, 2.0, 0.03, 0.001

nouter, ninner = 100, 50
hclose, rclose, relax = 1e-6, 1e-6, 1.0


def get_model(ws, name, extras):
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

    # -- LAK: one lake, one horizontal connection, optionally one manning
    #    outlet and an auxiliary variable
    flopy.mf6.ModflowGwflak(
        gwf,
        auxiliary=["AUX1"] if extras else None,
        nlakes=1,
        noutlets=1 if extras else 0,
        packagedata=[[0, lake_strt, 1, aux_val] if extras else [0, lake_strt, 1]],
        connectiondata=[
            [0, 0, (0, 0, lake_col), "horizontal", 0.1, botm, top, delr / 2.0, delc]
        ],
        outlets=[[0, 0, -1, "manning", invert, width, rough, slope]]
        if extras
        else None,
        perioddata={
            0: [
                [0, "rainfall", rainfall],
                [0, "evaporation", evaporation],
                [0, "runoff", runoff],
                [0, "inflow", inflow],
                [0, "withdrawal", withdrawal],
            ],
        },
        pname=lak_pname,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return sim


def build_models(idx, test):
    name = cases[idx]
    sim = get_model(test.workspace, name, with_extras[idx])
    # comparison sim in a "libmf6" subdir: the framework runs this one via the API
    sim_api = get_model(os.path.join(test.workspace, "libmf6"), name, with_extras[idx])
    return sim, sim_api


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    name = cases[idx].upper()
    extras = with_extras[idx]
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
        # advance one step so period (RAINFALL/WITHDRAWAL) data is populated
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

        pkg_prefix = f"{name}/{lak_pname}"
        input_prefix = f"__INPUT__/{name}/{lak_pname}"
        checked = [r for r in renames if extras or r[0] not in outlet_names]

        # -- every old name resolves and shares memory with the current name
        values = {}
        for old_name, new_name, in_input in checked:
            old_val = mf6.get_value(mf6.get_var_address(old_name, name, lak_pname))
            prefix = input_prefix if in_input else pkg_prefix
            new_val = mf6.get_value(f"{prefix}/{new_name}")
            assert np.array_equal(old_val, new_val), (
                f"alias '{old_name}' and '{new_name}' differ: {old_val} vs {new_val}"
            )
            values[old_name] = np.ravel(old_val)

        # -- old names carry the model input
        expected = {
            "RAINFALL": rainfall,
            "EVAPORATION": evaporation,
            "RUNOFF": runoff,
            "INFLOW": inflow,
            "WITHDRAWAL": withdrawal,
            "CELLID": lake_col + 1,  # reduced node number of the connected cell
            "IFORCEFB": 0,
        }
        if extras:
            expected.update(
                {
                    "OUTINVERT": invert,
                    "OUTWIDTH": width,
                    "OUTROUGH": rough,
                    "OUTSLOPE": slope,
                    "LAUXVAR": aux_val,
                }
            )
        for old_name, value in expected.items():
            assert np.allclose(values[old_name], value), (
                f"old name '{old_name}' value {values[old_name]} != {value}"
            )

        if not extras:
            # -- as in 6.8.1: LAUXVAR is empty and outlet names do not exist
            assert values["LAUXVAR"].size == 0, (
                f"old name 'LAUXVAR' is not empty: {values['LAUXVAR']}"
            )
            for old_name in outlet_names:
                tag = mf6.get_var_address(old_name, name, lak_pname)
                try:
                    mf6.get_value(tag)
                except Exception:
                    continue
                raise AssertionError(f"'{old_name}' resolves without outlets")

        # -- prove shared storage: write through the old name, read back the
        #    input-context variable
        tag = mf6.get_var_address("WITHDRAWAL", name, lak_pname)
        mf6.set_value(tag, np.array([withdrawal + 1.0]))
        after = mf6.get_value(f"{input_prefix}/WITHDRAWAL")
        assert np.isclose(after[0], withdrawal + 1.0), (
            "writing through deprecated name 'WITHDRAWAL' was not reflected in "
            f"input-context 'WITHDRAWAL' (got {after[0]})"
        )
        mf6.set_value(tag, np.array([withdrawal]))

    except Exception as e:
        print("Deprecated-name resolution check failed:")
        print(f"{type(e).__name__}: {e}")
        return False, open(output_file_path).readlines()

    try:
        mf6.finalize()
    except Exception:
        return False, open(output_file_path).readlines()

    # -- each old name used above is reported once; absent names are not
    expected = {f"{old}@{name}/{lak_pname}": 1 for old, _, _ in checked}
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

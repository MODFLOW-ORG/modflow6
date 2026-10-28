"""
Verify that deprecated advanced transport variable names remain resolvable
through the API.

IDM integration moved the advanced transport package (APT) feature and period
arrays from the package memory path into the input context, renaming them. The
APT packages register aliases at the package path so external programs
(modflowapi, coupled models) can still address these variables by their old
names. The models are reused from existing tests:

    case  model test                  package  package-specific old names
    sft   test_gwt_sft01              SFT      CONCRAIN, CONCEVAP, CONCROFF, CONCIFLW
    mwt   test_gwt_mwt01              MWT      CONCRATE
    uzt   test_gwt_uzt01              UZT      CONCINFL, CONCUZET
    lke   test_gwe_lke_conduction     LKE      TEMPRAIN, TEMPEVAP, TEMPROFF, TEMPIFLW
    lkt   test_gwt_lkt01              LKT      CONCRAIN, CONCEVAP, CONCROFF, CONCIFLW
    sfe   test_gwe_sfe01sto           SFE      TEMPRAIN, TEMPEVAP, TEMPROFF, TEMPIFLW
    mwe   test_gwe_mwe_conn_status    MWE      TEMPRATE
    uze   test_gwe_uze00              UZE      TEMPINFL, TEMPUZET

All packages also alias CONCFEAT to CONCENTRATION (GWT) or TEMPERATURE (GWE)
and LAUXVAR to AUX. SFT, MWT and LKT have auxiliary variables; the others do
not, so LAUXVAR is empty there, as in 6.8.1.

For each: the old name resolves and shares memory with the current
input-context variable, and each old name used is reported once as deprecated
in the simulation listing file.
"""

import os
import re
from importlib import import_module
from types import SimpleNamespace

import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

# -- case: (model builder module, model name, package name, auxiliary values,
#           [(old name, current name), ...])
apt_cases = {
    "libapt_sft": (
        "test_gwt_sft01",
        "GWT_SFT_01",
        "SFR-1",
        [99.0, 999.0],
        [
            ("CONCFEAT", "CONCENTRATION"),
            ("LAUXVAR", "AUX"),
            ("CONCRAIN", "RAINFALL"),
            ("CONCEVAP", "EVAPORATION"),
            ("CONCROFF", "RUNOFF"),
            ("CONCIFLW", "INFLOW"),
        ],
    ),
    "libapt_mwt": (
        "test_gwt_mwt01",
        "GWT_MWT_01",
        "MAW-1",
        [99.0, 999.0],
        [
            ("CONCFEAT", "CONCENTRATION"),
            ("LAUXVAR", "AUX"),
            ("CONCRATE", "RATE"),
        ],
    ),
    "libapt_uzt": (
        "test_gwt_uzt01",
        "GWT_UZT01A",
        "UZT-1",
        [],
        [
            ("CONCFEAT", "CONCENTRATION"),
            ("LAUXVAR", "AUX"),
            ("CONCINFL", "INFILTRATION"),
            ("CONCUZET", "UZET"),
        ],
    ),
    "libapt_lke": (
        "test_gwe_lke_conduction",
        "GWE-LKE-CONDUCTN",
        "LKE-1",
        [],
        [
            ("CONCFEAT", "TEMPERATURE"),
            ("LAUXVAR", "AUX"),
            ("TEMPRAIN", "RAINFALL"),
            ("TEMPEVAP", "EVAPORATION"),
            ("TEMPROFF", "RUNOFF"),
            ("TEMPIFLW", "EXT_INFLOW"),
        ],
    ),
    "libapt_lkt": (
        "test_gwt_lkt01",
        "GWT_LKT_01",
        "LKT-1",
        [99.0, 999.0],
        [
            ("CONCFEAT", "CONCENTRATION"),
            ("LAUXVAR", "AUX"),
            ("CONCRAIN", "RAINFALL"),
            ("CONCEVAP", "EVAPORATION"),
            ("CONCROFF", "RUNOFF"),
            ("CONCIFLW", "EXT_INFLOW"),
        ],
    ),
    "libapt_sfe": (
        "test_gwe_sfe01sto",
        "SFE_01STO_GWE",
        "SFR-1",
        [],
        [
            ("CONCFEAT", "TEMPERATURE"),
            ("LAUXVAR", "AUX"),
            ("TEMPRAIN", "RAINFALL"),
            ("TEMPEVAP", "EVAPORATION"),
            ("TEMPROFF", "RUNOFF"),
            ("TEMPIFLW", "INFLOW"),
        ],
    ),
    "libapt_mwe": (
        "test_gwe_mwe_conn_status",
        "GWE-MWE_CS",
        "MWE-1",
        [],
        [
            ("CONCFEAT", "TEMPERATURE"),
            ("LAUXVAR", "AUX"),
            ("TEMPRATE", "RATE"),
        ],
    ),
    "libapt_uze": (
        "test_gwe_uze00",
        "GWE_UZE00",
        "UZE-1",
        [],
        [
            ("CONCFEAT", "TEMPERATURE"),
            ("LAUXVAR", "AUX"),
            ("TEMPINFL", "INFILTRATION"),
            ("TEMPUZET", "UZET"),
        ],
    ),
}
cases = list(apt_cases)


def build_models(idx, test):
    module = import_module(apt_cases[cases[idx]][0])
    sim, _ = module.build_models(0, test)
    # -- copy in a "libmf6" subdir, run through the API by check_output
    api_ws = SimpleNamespace(workspace=test.workspace / "libmf6")
    sim_api, _ = module.build_models(0, api_ws)
    return sim, sim_api


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    _, model, pname, aux_vals, renames = apt_cases[cases[idx]]
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
        # -- advance one step so period data is populated
        mf6.update()

        input_prefix = f"__INPUT__/{model}/{pname}"

        # -- every old name resolves and shares memory with the current name
        for old_name, new_name in renames:
            old_val = mf6.get_value(mf6.get_var_address(old_name, model, pname))
            new_val = mf6.get_value(f"{input_prefix}/{new_name}")
            assert np.array_equal(old_val, new_val), (
                f"alias '{old_name}' and '{new_name}' differ: {old_val} vs {new_val}"
            )

        # -- auxiliary values, empty without auxiliary variables as in 6.8.1;
        #    the API returns one row per feature
        lauxvar = mf6.get_value(mf6.get_var_address("LAUXVAR", model, pname))
        if aux_vals:
            assert np.allclose(lauxvar[0], aux_vals), (
                f"old name 'LAUXVAR' feature-1 values {lauxvar[0]} != {aux_vals}"
            )
        else:
            assert lauxvar.size == 0, f"old name 'LAUXVAR' is not empty: {lauxvar}"

        # -- prove shared storage: write through the old name, read back the
        #    input-context variable
        _, feat_name = renames[0]
        tag = mf6.get_var_address("CONCFEAT", model, pname)
        concfeat = mf6.get_value(tag)
        mf6.set_value(tag, concfeat + 1.0)
        after = mf6.get_value(f"{input_prefix}/{feat_name}")
        assert np.allclose(after, concfeat + 1.0), (
            f"writing through deprecated name 'CONCFEAT' was not reflected in "
            f"input-context '{feat_name}'"
        )
        mf6.set_value(tag, concfeat)

    except Exception as e:
        print("Deprecated-name resolution check failed:")
        print(f"{type(e).__name__}: {e}")
        return False, open(output_file_path).readlines()

    try:
        mf6.finalize()
    except Exception:
        return False, open(output_file_path).readlines()

    # -- each old name used above is reported once
    expected = {f"{old}@{model}/{pname}": 1 for old, _ in renames}
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


def check_output(idx, test):
    # -- run the "libmf6" copy through the API; output comparison with the
    #    executable run is not the purpose of this test
    success, buff = api_func(test.targets["libmf6"], idx, test.workspace / "libmf6")
    assert success, "".join(buff[-20:])


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
    )
    test.run()

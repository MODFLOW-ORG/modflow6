import os
from pathlib import Path
from shutil import copytree

import flopy
import pytest
from compare import (
    Comparison,
    detect_comparison,
    setup_comparison,
    setup_mf5to6,
    setup_simulation,
)
from framework import TestFramework
from modflow_devtools.models import DEFAULT_REGISTRY, LocalRegistry

PREFIXES = ["mf6/test", "mf6/large", "mf2005"]
EXCLUDE = {
    "alt_model",
    "test205_gwtbuy-henrytidal",
    # todo reinstate after 6.5.0 release
    "test001d_Tnewton",
    # remove tests with nwt usg conductance weighting
    "test006_gwf3_gnc_nr_dev",
    "test006_gwf3_nr_dev",
    "test014_NWTP3High_dev",
    "test015_KeatingLike_disu_dev",
    "test041_flowdivert_nr_dev",
    "test016_Keating_disu_dev",
    "test053_npf-a-nwt_dev",
    "test053_npf-b-nwt_dev",
    # todo reinstate after resolving convergence failure
    "test014_NWTP3Low_dev",
    "test1002_biscqtg_disv_gnc_nr_dev",
    "test1002_biscqtg_disv_nr_MD_dev",
    "test1002_biscqtg_disv_nr_RCM_dev",
    "test1002_biscqtg_disv_nr_dev",
    "testWetDry/mf2005",
}


def pytest_generate_tests(metafunc):
    if "model_name" in metafunc.fixturenames:
        models_paths = metafunc.config.getoption("--models-path")
        models_paths = [
            Path(p).expanduser().resolve().absolute() for p in models_paths or []
        ]
        registry = LocalRegistry() if any(models_paths) else DEFAULT_REGISTRY
        registry_type = type(registry).__name__.lower().replace("registry", "")
        metafunc.parametrize("registry", [registry], ids=[registry_type])
        models = []
        if "local" in registry_type:
            namefile_pattern = (
                metafunc.config.getoption("--namefile-pattern") or "mfsim.nam"
            )
            for path in models_paths:
                registry.index(path, namefile=namefile_pattern)
            models.extend(registry.models.keys())
        else:
            for model_prefix in PREFIXES:
                models.extend(
                    [m for m in registry.models.keys() if m.startswith(model_prefix)]
                )
        models = sorted(models)
        metafunc.parametrize("model_name", models, ids=models)


@pytest.mark.slow
@pytest.mark.external
@pytest.mark.regression
def test_model(
    registry,
    model_name,
    tmp_path,
    markers,
    targets,
    function_tmpdir,
    original_regression,
):
    exclude = any(s in model_name for s in EXCLUDE)
    devonly = "dev" in model_name and "not developmode" in markers
    large = "large" in model_name and "not large" in markers
    if exclude or devonly or large:
        reason = "excluded" if exclude else "developmode only"
        pytest.skip(f"Skipping: {model_name} ({reason})")

    # TODO: avoid this intermediate copy? it's needed
    # because the simulation workspace should only be
    # a subset of all the model directory contents in
    # some cases. maybe allow filtering in `copy_to`?
    registry.copy_to(tmp_path, model_name)
    if "mfsim.nam" in os.listdir(tmp_path):
        src_workspace = tmp_path
        mf6_workspace = function_tmpdir
    else:
        # mf2005 model, run the mf5to6 converter
        mf6_workspace = function_tmpdir / "mf6"
        src_workspace = function_tmpdir / "mf5to6"
        npth = setup_mf5to6(tmp_path, src_workspace)
        nam = os.path.basename(npth)
        exe = os.path.abspath(targets["mf5to6"])
        print("MODFLOW 5 to 6 converter run for", nam, "using executable", exe)
        success, _ = flopy.run_model(
            exe,
            nam,
            model_ws=src_workspace,
            normal_msg="Program terminated normally",
            cargs="mf6",
        )
        assert success

    setup_simulation(src=src_workspace, dst=mf6_workspace)
    if (
        comparison := detect_comparison(tmp_path)
        if original_regression
        else Comparison.MF6_REGRESSION
    ) == Comparison.MF6_REGRESSION:
        copytree(src_workspace, mf6_workspace / comparison.value)
    else:
        setup_comparison(
            mf6_workspace,
            mf6_workspace / comparison.value,
            comparison.value,
            overwrite=True,
        )

    test = TestFramework(
        name=model_name,
        workspace=mf6_workspace,
        targets=targets,
        compare=comparison,
        verbose=False,
    )
    test.run()

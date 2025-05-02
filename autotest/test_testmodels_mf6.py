from pathlib import Path
from shutil import copytree

import pytest
from compare import (
    Comparison,
    detect_comparison,
    setup_comparison,
    setup_simulation,
)
from framework import TestFramework
from modflow_devtools.models import DEFAULT_REGISTRY, LocalRegistry

SKIP = [
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
]


def pytest_generate_tests(metafunc):
    paths = [
        Path(p).expanduser().resolve().absolute()
        for p in metafunc.config.getoption("--models-path") or []
    ]
    if any(paths):
        registry = LocalRegistry()
        for path in paths:
            registry.index(path)
        models = list(registry.models.keys())
        metafunc.parametrize("registry", [registry], ids=["local-registry"])
    else:
        registry = DEFAULT_REGISTRY
        models = [m for m in registry.models.keys() if m.startswith("mf6/test/")]
        metafunc.parametrize("registry", [registry], ids=["default-registry"])
    if "model_name" in metafunc.fixturenames:
        metafunc.parametrize("model_name", models, ids=models)


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
    skip = any(s in model_name for s in SKIP)
    devonly = "dev" in model_name and "not developmode" in markers
    if skip or devonly:
        reason = "excluded" if skip else "developmode only"
        pytest.skip(f"Skipping: {model_name} ({reason})")

    # TODO: avoid this intermediate copy? it's needed
    # because the simulation workspace should only be
    # a subset of all the model directory contents in
    # some cases. maybe allow filtering in `copy_to`?
    registry.copy_to(tmp_path, model_name)
    setup_simulation(src=tmp_path, dst=function_tmpdir)
    if (
        comparison := detect_comparison(tmp_path)
        if original_regression
        else Comparison.MF6_REGRESSION
    ) == Comparison.MF6_REGRESSION:
        copytree(function_tmpdir, function_tmpdir / comparison.value)
    else:
        setup_comparison(
            function_tmpdir,
            function_tmpdir / comparison.value,
            comparison.value,
            overwrite=True,
        )

    test = TestFramework(
        name=model_name,
        workspace=function_tmpdir,
        targets=targets,
        compare=comparison,
        verbose=False,
    )
    test.run()

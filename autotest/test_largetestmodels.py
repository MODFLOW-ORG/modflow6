from shutil import copytree

import modflow_devtools.models as models
import pytest
from compare import (
    Comparison,
    detect_comparison,
    setup_comparison,
    setup_simulation,
)
from framework import TestFramework

MODELS = [m for m in models.get_models().keys() if m.startswith("large/")]
SKIP = [
    "test1002_biscqtg_disv_gnc_nr_dev",
    "test1002_biscqtg_disv_nr_MD_dev",
    "test1002_biscqtg_disv_nr_RCM_dev",
    "test1002_biscqtg_disv_nr_dev",
]


@pytest.mark.regression
@pytest.mark.slow
@pytest.mark.parametrize("model_name", MODELS)
def test_model(
    model_name,
    tmp_path,
    function_tmpdir,
    markers,
    original_regression,
    targets,
):
    models.copy_to(tmp_path, model_name)

    # skip this one?
    skip = any(s in model_name for s in SKIP)
    devonly = "dev" in model_name and "not developmode" in markers
    if skip or devonly:
        reason = "excluded" if skip else "developmode only"
        pytest.skip(f"Skipping: {model_name} ({reason})")

    # setup test workspace and framework
    setup_simulation(src=tmp_path, dst=function_tmpdir)

    # setup comparison workspace
    if (
        compare := detect_comparison(tmp_path)
        if original_regression
        else Comparison.MF6_REGRESSION
    ) == Comparison.MF6_REGRESSION:
        copytree(function_tmpdir, function_tmpdir / compare.value)
    else:
        setup_comparison(
            function_tmpdir,
            function_tmpdir / compare.value,
            compare.value,
            overwrite=True,
        )

    # run the test
    test = TestFramework(
        name=model_name,
        workspace=function_tmpdir,
        targets=targets,
        compare=compare,
        verbose=False,
    )
    test.run()

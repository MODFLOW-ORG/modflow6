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

MODELS = [m for m in models.get_models().keys() if m.startswith("test/")]
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


@pytest.mark.regression
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

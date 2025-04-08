import os
from pathlib import Path
from shutil import copytree

import flopy
import modflow_devtools.models as models
import pytest
from compare import (
    Comparison,
    detect_comparison,
    get_namefiles,
    setup_comparison,
    setup_model,
    setup_simulation,
)
from framework import TestFramework

MODELS = [m for m in models.get_models().keys() if m.startswith("mf2005/")]
SKIP = ["alt_model", "Dry/mf2005"]


def setup_mf5to6(src, dst) -> Path:
    Path(dst).mkdir(exist_ok=True)
    lgrpth = None

    # determine if compare directory exists in directory or if mflgr control
    # file is in directory
    listdir = os.listdir(src)
    for value in listdir:
        fpth = os.path.join(src, value)
        if os.path.isfile(fpth):
            ext = os.path.splitext(fpth)[1]
            if ".lgr" in ext.lower():
                lgrpth = fpth

    print(f"Copying files to target workspace: {dst}")
    # copy lgr files to working directory
    if lgrpth is not None:
        npth = lgrpth
        setup_model(lgrpth, dst)
    # copy MODFLOW-2005, MODFLOW-NWT, or MODFLOW-USG files to working directory
    else:
        npths = get_namefiles(src)
        if len(npths) < 1:
            msg = f"No name files in source workspace: {src}"
            print(msg)
            assert False
        npth = npths[0]
        setup_model(npth, dst)

    return Path(npth)


@pytest.mark.external
@pytest.mark.regression
@pytest.mark.parametrize("model_name", MODELS)
def test_model(
    model_name,
    tmp_path,
    function_tmpdir,
    original_regression,
    targets,
):
    models.copy_to(tmp_path, model_name)

    if any(s in model_name for s in SKIP):
        pytest.skip(f"Skipping: {model_name} (excluded)")

    # run the mf5to6 converter
    mf5to6_workspace = function_tmpdir / "mf5to6"
    npth = setup_mf5to6(tmp_path, mf5to6_workspace)
    nam = os.path.basename(npth)
    exe = os.path.abspath(targets["mf5to6"])
    print("MODFLOW 5 to 6 converter run for", nam, "using executable", exe)
    success, _ = flopy.run_model(
        exe,
        nam,
        model_ws=mf5to6_workspace,
        normal_msg="Program terminated normally",
        cargs="mf6",
    )
    assert success

    # setup mf6 workspace
    mf6_workspace = function_tmpdir / "mf6"
    setup_simulation(src=mf5to6_workspace, dst=mf6_workspace)
    if (
        comparison := (
            detect_comparison(mf5to6_workspace)
            if original_regression
            else Comparison.MF6_REGRESSION
        )
    ) == Comparison.MF6_REGRESSION:
        copytree(mf5to6_workspace, mf6_workspace / comparison.value)
    else:
        setup_comparison(
            mf5to6_workspace, mf6_workspace, comparison.value, overwrite=True
        )

    # run the test
    test = TestFramework(
        name=model_name,
        workspace=mf6_workspace,
        targets=targets,
        compare=comparison,
        verbose=False,
    )
    test.run()

import os
from pathlib import Path
from shutil import copytree

import flopy
import pytest
from compare import (
    detect_comparison,
    get_namefiles,
    setup_comparison,
    setup_model,
    setup_simulation,
)
from framework import TestFramework

excluded_models = ["alt_model", "mf2005"]


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


@pytest.mark.repo
@pytest.mark.regression
def test_model(
    function_tmpdir,
    original_regression,
    targets,
    # https://modflow-devtools.readthedocs.io/en/latest/md/fixtures.html#modflow-2005-test-models
    test_model_mf5to6,
):
    model_path = test_model_mf5to6.parent
    model_name = model_path.name
    if model_name in excluded_models:
        pytest.skip(f"Skipping: {model_name} (excluded)")

    # run the mf5to6 converter
    mf5to6_workspace = function_tmpdir / "mf5to6"
    npth = setup_mf5to6(model_path, mf5to6_workspace)
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

    # setup mf6 workspace and framework
    mf6_workspace = function_tmpdir / "mf6"
    setup_simulation(src=mf5to6_workspace, dst=mf6_workspace)
    compare = (
        detect_comparison(mf5to6_workspace) if original_regression else "mf6_regression"
    )
    test = TestFramework(
        name=model_path.name,
        workspace=mf6_workspace,
        targets=targets,
        compare=compare,
        verbose=False,
    )

    # setup comparison workspace
    if compare == "mf6_regression":
        copytree(mf5to6_workspace, mf6_workspace / compare)
    else:
        setup_comparison(mf5to6_workspace, mf6_workspace, compare, overwrite=True)

    # run the test
    test.run()

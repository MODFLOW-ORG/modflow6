"""
Test to make sure that mf6 is failing with the correct error messages.  This
test script is set up to be extensible so that simple models can be created
very easily and tested with different options to succeed or fail correctly.
"""

import subprocess

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs


def run_mf6(argv, ws):
    buff = []
    proc = subprocess.Popen(
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, error = proc.communicate()
    if result is not None:
        c = result.decode("utf-8")
        c = c.rstrip("\r\n")
        print(f"{c}")
        buff.append(c)

    return proc.returncode, buff


def run_mf6_error(ws, exe, err_str_list):
    returncode, buff = run_mf6([exe], ws)
    msg = "mf terminated with error"
    if returncode != 0:
        if not isinstance(err_str_list, list):
            err_str_list = [err_str_list]
        for err_str in err_str_list:
            err = any(err_str in s for s in buff)
            if err:
                raise RuntimeError(msg)
            else:
                msg += " but did not print correct error message."
                msg += f'  Correct message should have been "{err_str}"'
                raise ValueError(msg)


def get_minimal_gwf_simulation(
    ws,
    exe,
    name="test",
    simkwargs=None,
    simnamefilekwargs=None,
    tdiskwargs=None,
    gwfkwargs=None,
    imskwargs=None,
    diskwargs=None,
    disukwargs=None,
    ickwargs=None,
    npfkwargs=None,
    chdkwargs=None,
):
    if simkwargs is None:
        simkwargs = {}
    if tdiskwargs is None:
        tdiskwargs = {}
    if gwfkwargs is None:
        gwfkwargs = {}
        gwfkwargs["modelname"] = name
    if imskwargs is None:
        imskwargs = {"print_option": "SUMMARY"}
    if diskwargs is None and disukwargs is None:
        diskwargs = {}
        diskwargs["nlay"] = 5
        diskwargs["nrow"] = 5
        diskwargs["ncol"] = 5
        diskwargs["top"] = 0
        diskwargs["botm"] = [-1, -2, -3, -4, -5]
    if ickwargs is None:
        ickwargs = {}
    if npfkwargs is None:
        npfkwargs = {}
    if chdkwargs is None:
        chdkwargs = {}
        nl = diskwargs["nlay"]
        nr = diskwargs["nrow"]
        nc = diskwargs["ncol"]
        chdkwargs["stress_period_data"] = {
            0: [[(0, 0, 0), 0], [(0, nr - 1, nc - 1), 1]]
        }
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name=exe, sim_ws=ws, **simkwargs
    )
    if simnamefilekwargs is not None:
        for k in simnamefilekwargs:
            sim.name_file.__setattr__(k, simnamefilekwargs[k])
    tdis = flopy.mf6.ModflowTdis(sim, **tdiskwargs)
    gwf = flopy.mf6.ModflowGwf(sim, **gwfkwargs)
    ims = flopy.mf6.ModflowIms(sim, **imskwargs)
    if diskwargs is not None:
        dis = flopy.mf6.ModflowGwfdis(gwf, **diskwargs)
    elif disukwargs is not None:
        disu = flopy.mf6.ModflowGwfdisu(gwf, **disukwargs)
    ic = flopy.mf6.ModflowGwfic(gwf, **ickwargs)
    npf = flopy.mf6.ModflowGwfnpf(gwf, **npfkwargs)
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, **chdkwargs)
    return sim


def test_simple_model_success(function_tmpdir, targets):
    mf6 = targets["mf6"]

    # test a simple model to make sure it runs and terminates correctly
    sim = get_minimal_gwf_simulation(str(function_tmpdir), mf6)
    sim.write_simulation()
    returncode, buff = run_mf6([mf6], str(function_tmpdir))
    assert returncode == 0, "mf6 failed for simple model."

    final_message = "Normal termination of simulation."
    failure_message = f'mf6 did not terminate with "{final_message}"'
    assert final_message in buff[-1], failure_message


def test_input_bound(function_tmpdir, targets):
    mf6 = targets["mf6"]

    # test a simple model to make sure it runs and terminates correctly
    sim = get_minimal_gwf_simulation(str(function_tmpdir), mf6)
    sim.write_simulation()

    # update chd to underspecify maxbound
    with open(function_tmpdir / "test.chd", "w") as f:
        f.write("BEGIN options\n")
        f.write("END options\n\n")
        f.write("BEGIN dimensions\n")
        f.write("  MAXBOUND  1\n")
        f.write("END dimensions\n\n")
        f.write("BEGIN period  1\n")
        f.write("  1 1 1 0\n")
        f.write("  1 5 5 1\n")
        f.write("END period  1\n")

    with pytest.raises(RuntimeError):
        # make sure error is set when input dimension is too small
        err_str = "Input error: line count exceeds input dimension. Expected rows=1."
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_empty_folder(function_tmpdir, targets):
    mf6 = targets["mf6"]
    with pytest.raises(RuntimeError):
        # make sure mf6 fails when there is no simulation name file
        err_str = "mfsim.nam is not present in working directory."
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_sim_errors(function_tmpdir, targets):
    mf6 = targets["mf6"]

    with pytest.raises(RuntimeError):
        # verify that the correct number of errors are reported
        chdkwargs = {}
        chdkwargs["stress_period_data"] = {0: [[(0, 0, 0), 0.0] for i in range(10)]}
        sim = get_minimal_gwf_simulation(
            str(function_tmpdir), exe=mf6, chdkwargs=chdkwargs
        )
        sim.write_simulation()
        err_str = ["1. Cell is already a constant head ((1,1,1))."]
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_sim_maxerrors(function_tmpdir, targets):
    mf6 = targets["mf6"]

    with pytest.raises(RuntimeError):
        # verify that the maxerrors keyword gives the correct error output
        simnamefilekwargs = {}
        simnamefilekwargs["maxerrors"] = 5
        chdkwargs = {}
        chdkwargs["stress_period_data"] = {0: [[(0, 0, 0), 0.0] for i in range(10)]}
        sim = get_minimal_gwf_simulation(
            str(function_tmpdir),
            exe=mf6,
            simnamefilekwargs=simnamefilekwargs,
            chdkwargs=chdkwargs,
        )
        sim.write_simulation()
        err_str = [
            "5. Cell is already a constant head ((1,1,1)).",
            "5 additional errors detected but not printed.",
            "UNIT ERROR REPORT:",
            "1. ERROR OCCURRED WHILE READING FILE 'test.chd'",
        ]
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_disu_errors(function_tmpdir, targets):
    mf6 = targets["mf6"]

    with pytest.raises(RuntimeError):
        disukwargs = get_disu_kwargs(3, 3, 3, np.ones(3), np.ones(3), 0, [-1, -2, -3])
        top = disukwargs["top"]
        bot = disukwargs["bot"]
        top[9] = 2.0
        bot[9] = 1.0
        sim = get_minimal_gwf_simulation(
            str(function_tmpdir),
            exe=mf6,
            disukwargs=disukwargs,
            chdkwargs={"stress_period_data": [[]]},
        )
        sim.write_simulation()
        err_str = [
            "1. Top elevation (    2.00000    ) for cell 10 is above bottom elevation (",  # noqa
            "-1.00000    ) for cell 1. Based on node numbering rules cell 10 must be",
            "below cell 1.",
            "UNIT ERROR REPORT:1. ERROR OCCURRED WHILE READING FILE './test.disu'",
        ]
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_solver_fail(function_tmpdir, targets):
    mf6 = targets["mf6"]

    with pytest.raises(RuntimeError):
        # test failed to converge
        imskwargs = {"inner_maximum": 1, "outer_maximum": 2}
        sim = get_minimal_gwf_simulation(
            str(function_tmpdir), exe=mf6, imskwargs=imskwargs
        )
        sim.write_simulation()
        err_str = [
            "Simulation convergence failure occurred 1 time(s).",
            "Premature termination of simulation.",
        ]
        run_mf6_error(str(function_tmpdir), mf6, err_str)


def test_fail_continue_success(function_tmpdir, targets):
    mf6 = targets["mf6"]

    # test continue but failed to converge
    tdiskwargs = {"nper": 1, "perioddata": [(10.0, 10, 1.0)]}
    imskwargs = {"inner_maximum": 1, "outer_maximum": 2}
    sim = get_minimal_gwf_simulation(
        str(function_tmpdir),
        exe=mf6,
        imskwargs=imskwargs,
        tdiskwargs=tdiskwargs,
    )
    sim.name_file.continue_ = True
    sim.write_simulation()
    returncode, buff = run_mf6([mf6], str(function_tmpdir))
    assert returncode == 0, "mf6 failed for simple model."

    final_message = "Simulation convergence failure occurred 10 time(s)."
    failure_message = f'mf6 did not terminate with "{final_message}"'
    assert final_message in buff[0], failure_message

    final_message = "Normal termination of simulation."
    failure_message = f'mf6 did not terminate with "{final_message}"'
    assert final_message in buff[0], failure_message

import os
import subprocess
import sys

import flopy
import pytest

name = "gwf_ret_codes01"
app = "mf6"
if sys.platform.lower() == "win32":
    app += ".exe"


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


def get_sim(ws, exe, idomain, continue_flag=False, nouter=500):
    # static model data
    # temporal discretization
    nper = 1
    tdis_rc = [(1.0, 1, 1.0)]

    # spatial discretization data
    nlay, nrow, ncol = 1, 10, 10
    delr, delc = 1000.0, 1000.0
    top = 0.0
    botm = [-100]
    strt = 0.0

    # calculate hk
    hk = 1.0

    # solver options
    ninner = 300
    hclose, rclose, relax = 1e-9, 1e-6, 1.0
    newtonoptions = "NEWTON"
    imsla = "BICGSTAB"

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=exe,
        sim_ws=ws,
        continue_=continue_flag,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions=newtonoptions, save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=1, k=hk)

    # constant head
    if idomain > 0:
        c6 = []
        ccol = [1, ncol - 1]
        for j in ccol:
            c6.append([(0, nrow - 1, j), strt])
        c6 = [[0, 0, 0, 1.0], [0, nrow - 1, ncol - 1, 0.0]]
        cd6 = {0: c6}
        maxchd = len(cd6[0])
        chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=cd6, maxbound=maxchd)

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    # write the input files
    sim.write_simulation()
    return sim


def normal_termination(dir, exe):
    ws = os.path.join(dir, "normal_termination")

    # get the simulation
    sim = get_sim(ws, exe, idomain=1)

    # write the input files
    sim.write_simulation()

    # run the simulation
    returncode, buff = run_mf6([exe], ws)
    if returncode != 0:
        msg = (
            "The run should have been successful but it terminated "
            "with non-zero returncode"
        )
        raise ValueError(msg)


def converge_fail_continue(dir, exe):
    ws = os.path.join(dir, "converge_fail_continue")

    # get the simulation
    sim = get_sim(ws, exe, idomain=1, continue_flag=True, nouter=1)

    # write the input files
    sim.write_simulation()

    # run the simulation
    returncode, buff = run_mf6([exe], ws)
    msg = (
        "The run should have been successful even though it failed, because"
        " the continue flag was set.  But a non-zero error code was "
        f"found: {returncode}"
    )
    assert returncode == 0, msg


def converge_fail_nocontinue(dir, exe):
    ws = os.path.join(dir, "converge_fail_nocontinue")

    with pytest.raises(RuntimeError):
        # get the simulation
        sim = get_sim(ws, exe, idomain=1, continue_flag=False, nouter=1)

        # write the input files
        sim.write_simulation()

        # run the simulation
        returncode, buff = run_mf6([exe], ws)
        msg = "This run should fail with a returncode of 1"
        if returncode == 1:
            raise RuntimeError(msg)


def idomain_runtime_error(dir, exe):
    ws = os.path.join(dir, "idomain_runtime_error")

    with pytest.raises(RuntimeError):
        # get the simulation
        sim = get_sim(ws, exe, idomain=0)

        # write the input files
        sim.write_simulation()

        # run the simulation
        returncode, buff = run_mf6([exe], ws)
        msg = f"could not run {sim.name}"
        if returncode != 0:
            err_str = "Ensure IDOMAIN array has some"
            err = any(err_str in s for s in buff)
            if err:
                raise RuntimeError(msg)
            else:
                msg += " but IDOMAIN ARRAY ERROR not returned"
                raise ValueError(msg)


def unknown_keyword_error(dir, exe):
    with pytest.raises((RuntimeError, ValueError)):
        returncode, buff = run_mf6([exe, "--unknown_keyword"], dir)
        msg = "could not run unknown_keyword"
        if returncode != 0:
            err_str = f"{app}: illegal option"
            err = any(err_str in s for s in buff)
            if err:
                raise RuntimeError(msg)
            else:
                msg += f" but {err_str} not returned"
                raise ValueError(msg)


def run_argv(arg, return_str, tempdir, exe):
    returncode, buff = run_mf6([exe, arg], tempdir)
    if returncode == 0:
        found_str = any(return_str in s for s in buff)
        if not found_str:
            msg = f"{arg} keyword did not return {return_str}"
            raise ValueError(msg)
    else:
        msg = f"could not run with command line argument {arg}"
        raise RuntimeError(msg)


def help_argv(dir, exe):
    for arg in ["-h", "--help", "-?"]:
        return_str = f"{app} [options]     retrieve program information"
        run_argv(arg, return_str, dir, exe)


def version_argv(dir, exe):
    for arg in ["-v", "--version"]:
        return_str = f"{app}: 6"
        run_argv(arg, return_str, dir, exe)


def develop_argv(dir, exe):
    for arg in ["-dev", "--develop"]:
        return_str = f"{app}: develop version"
        run_argv(arg, return_str, dir, exe)


def compiler_argv(dir, exe):
    for arg in ["-c", "--compiler"]:
        return_str = f"{app}: MODFLOW 6 compiled"
        run_argv(arg, return_str, dir, exe)


@pytest.mark.parametrize(
    "fn",
    (
        "idomain_runtime_error",
        "unknown_keyword_error",
        "normal_termination",
        "converge_fail_nocontinue",
        "help_argv",
        "version_argv",
        "develop_argv",
        "compiler_argv",
    ),
)
def test_main(fn, function_tmpdir, targets):
    mf6 = targets["mf6"]
    eval(fn)(function_tmpdir, mf6)

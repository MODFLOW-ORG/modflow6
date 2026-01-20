"""
Test the bmi which is used to calculate a recharge rate that results in a
simulated head in the center of the model domain to be equal to the
simulated head in the non-bmi simulation.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework
from modflow_devtools.markers import requires_pkg

cases = ["libgwf_rch02"]

# recharge package name
rch_pname = "RCH-1"


# average recharge rate
avg_rch = 0.001
drch = 1e-6 * avg_rch

# calculate recharge rates
dx = 1 / 20
rad = np.arange(0, 1 + dx, dx) * 2.0 * np.pi
f = np.sin(rad)
rch_rates = avg_rch + f * avg_rch

# temporal discretization
nper = rch_rates.shape[0]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((1.0, 1, 1))

# model spatial dimensions
nlay, nrow, ncol = 1, 11, 11

# cell spacing
delr = 10.0
delc = 10.0
area = delr * delc

# top of the aquifer
top = 10.0

# bottom of the aquifer
botm = 0.0

# hydraulic conductivity
hk = 1.0

# starting head
strt = 5.0

# build chd stress period data
chd_spd = {0: [[(0, 0, 0), strt], [(0, nrow - 1, ncol - 1), strt]]}

# build recharge spd
rch_spd = {}
for n in range(nper):
    rch_spd[n] = rch_rates[n]

# solver data
nouter, ninner = 100, 100
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def get_model(ws, name, exe, rech=rch_spd):
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="all",
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="SIMPLE",
        under_relaxation_gamma=0.98,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )

    # create gwf model
    newtonoptions = "NEWTON UNDER_RELAXATION"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        newtonoptions=newtonoptions,
        modelname=name,
        print_input=True,
        save_flows=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=hk)

    # chd file
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)

    # recharge file
    rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=rech, pname=rch_pname)

    # gwf observations
    onam = f"{name}.head.obs"
    cnam = onam + ".csv"
    obs_recarray = {cnam: [("h1_6_6", "HEAD", (0, 5, 5))]}
    gwfobs = flopy.mf6.ModflowUtlobs(
        gwf,
        print_input=True,
        filename=onam,
        digits=20,
        continuous=obs_recarray,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name, "mf6")

    # build comparison model
    ws = os.path.join(test.workspace, "libmf6")
    mc = get_model(ws, name, "mf6", rech=0.0)

    return sim, mc


def run_perturbation(mf6, max_iter, recharge, tag, rch):
    mf6.prepare_solve()
    kiter = 0
    while kiter < max_iter:
        # update recharge
        recharge[:] = rch
        mf6.set_value(tag, recharge)
        # solve with updated well rate
        has_converged = mf6.solve()
        kiter += 1
        if has_converged:
            break
    return has_converged


def api_func(exe, idx, model_ws=None):
    from modflowapi import ModflowApi

    print("\nBMI implementation test:")

    name = cases[idx].upper()
    init_wd = os.path.abspath(os.getcwd())
    if model_ws is not None:
        os.chdir(model_ws)

    output_file_path = os.path.join(model_ws, "mfsim.stdout")

    # get the observations from the standard run
    fpth = os.path.join("..", f"{cases[idx]}.head.obs.csv")
    hobs = np.genfromtxt(fpth, delimiter=",", names=True)["H1_6_6"]

    try:
        mf6 = ModflowApi(exe)
    except Exception as e:
        print("Failed to load " + str(exe))
        print("with message: " + str(e))
        return False, open(output_file_path).readlines()

    # initialize the model
    try:
        mf6.initialize()
    except:
        return False, open(output_file_path).readlines()

    # time loop
    current_time = mf6.get_current_time()
    end_time = mf6.get_end_time()

    # get pointer to simulated heads
    head_tag = mf6.get_var_address("X", name)
    head = mf6.get_value_ptr(head_tag)

    # maximum outer iterations
    mxit_tag = mf6.get_var_address("MXITER", "SLN_1")
    max_iter = mf6.get_value(mxit_tag)

    # get copy of recharge array
    rch_tag = mf6.get_var_address("RECHARGE", name, rch_pname)
    new_recharge = mf6.get_value(rch_tag).copy()

    # determine initial recharge value
    np.random.seed(0)
    rch = np.random.normal(1) * avg_rch

    # model time loop
    idx = 0
    while current_time < end_time:
        # target head
        htarget = hobs[idx]

        # get dt and prepare for non-linear iterations
        dt = mf6.get_time_step()
        mf6.prepare_time_step(dt)

        est_iter = 0
        while est_iter < 100:
            # base simulation loop
            has_converged = run_perturbation(mf6, max_iter, new_recharge, rch_tag, rch)
            if not has_converged:
                return False, open(output_file_path).readlines()
            h0 = head.reshape((nrow, ncol))[5, 5]
            r0 = h0 - htarget

            # perturbation simulation loop
            has_converged = run_perturbation(
                mf6, max_iter, new_recharge, rch_tag, rch + drch
            )
            if not has_converged:
                return False, open(output_file_path).readlines()
            h1 = head.reshape((nrow, ncol))[5, 5]
            r1 = h1 - htarget

            # calculate update terms
            dqdr = drch / (r0 - r1)
            dr = r1 * dqdr

            # evaluate if the estimation iterations need to continue
            if abs(r0) < 1e-5:
                msg = (
                    f"Estimation for time {current_time:5.1f}"
                    + f" converged in {est_iter:3d}"
                    + " iterations"
                    + f" -- final recharge={rch:10.5f}"
                    + f" residual={rch - rch_rates[idx]:10.2g}"
                )
                print(msg)
                break
            else:
                est_iter += 1
                rch += dr

        # solution with final estimated recharge for the timestep
        has_converged = run_perturbation(mf6, max_iter, new_recharge, rch_tag, rch)
        if not has_converged:
            return False, open(output_file_path).readlines()

        # finalize time step
        mf6.finalize_solve()

        # finalize time step and update time
        mf6.finalize_time_step()
        current_time = mf6.get_current_time()

        # increment counter
        idx += 1
    # cleanup
    try:
        mf6.finalize()
    except:
        return False, open(output_file_path).readlines()

    if model_ws is not None:
        os.chdir(init_wd)

    # cleanup and return
    return True, open(output_file_path).readlines()


@requires_pkg("modflowapi")
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
        api_func=lambda exe, ws: api_func(exe, idx, ws),
    )
    test.run()

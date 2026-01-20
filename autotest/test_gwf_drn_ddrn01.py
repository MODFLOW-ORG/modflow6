import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["drn_ddrn01a", "drn_ddrn01b"]
paktest = "drn"
budtol = 1e-2
ddir = "data"
newton = [False, True]

# static model data
# spatial discretization
nlay, nrow, ncol = 1, 1, 100
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
xlen = 1000.0
delc, delr = 1.0, xlen / ncol
top, botm = 10.0, 0.0

# temporal discretization
nper = 1
nyears = 10.0
period_data = [(10.0 * nyears, 100, 1.1)]

kh, sy, ss = 10.0, 0.1, 1e-5
h0, h1 = 9.0, 1.0
obs_recarray = {
    "head_obs.csv": [
        ("h1_1_1", "HEAD", (0, 0, 0)),
        ("h1_1_100", "HEAD", (0, 0, ncol - 1)),
    ]
}
# drain data
delev, ddrn = 0.0, h1
dcond = kh * delc * ddrn / (0.5 * delr)
drn_spd = {0: [[(0, 0, ncol - 1), delev, dcond, ddrn]]}
drn_obs = {"drn_obs.csv": [("d1_1_100", "DRN", (0, 0, ncol - 1))]}


def initial_conditions():
    x = np.arange(0, xlen - delr / 2, delr)
    return np.sqrt(h0**2 + x * (h1**2 - h0**2) / (xlen - delr))


def get_model(idx, ws, name):
    strt = initial_conditions()
    hdsfile = f"{name}.hds"
    if newton[idx]:
        newtonoptions = "NEWTON"
    else:
        newtonoptions = None

    # build the model
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=period_data)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        linear_acceleration="BICGSTAB",
        outer_maximum=200,
        inner_maximum=200,
        outer_dvclose=1e-6,
        inner_dvclose=1e-9,
        rcloserecord=[0.01, "strict"],
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions=newtonoptions)
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, k=kh, icelltype=1)
    sto = flopy.mf6.ModflowGwfsto(gwf, sy=sy, ss=ss, transient={0: True}, iconvert=1)
    drn = flopy.mf6.ModflowGwfdrn(
        gwf,
        auxiliary=["ddrn"],
        auxdepthname="ddrn",
        stress_period_data=drn_spd,
        print_input=True,
    )
    drn.obs.initialize(
        filename=f"{name}.drn.obs",
        digits=20,
        print_input=True,
        continuous=drn_obs,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
    obs_package = flopy.mf6.ModflowUtlobs(
        gwf, digits=20, print_input=True, continuous=obs_recarray
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=hdsfile,
        saverecord=[("HEAD", "ALL")],
        printrecord=[("BUDGET", "ALL")],
    )

    return sim


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = get_model(idx, ws, name)

    return sim, None


def check_output(idx, test):
    # MODFLOW 6 drain discharge results
    fpth = os.path.join(test.workspace, "drn_obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # MODFLOW 6 head results
    fpth = os.path.join(test.workspace, "head_obs.csv")
    try:
        th0 = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate the drain flux analytically
    xdiff = th0["H1_1_100"] - delev
    f = drain_smoothing(xdiff, ddrn, newton=newton[idx])
    tc0 = f * dcond * (delev - th0["H1_1_100"])

    # calculate maximum absolute error
    diff = tc["D1_1_100"] - tc0
    diffmax = np.abs(diff).max()
    dtol = 1e-6
    msg = f"maximum absolute discharge difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.disc.cmp.out")
    with open(fpth, "w") as f:
        line = f"{'TOTIM':>15s}"
        line += f" {'DRN':>15s}"
        line += f" {'UZF':>15s}"
        line += f" {'DIFF':>15s}"
        f.write(line + "\n")
        for i in range(diff.shape[0]):
            line = f"{tc['time'][i]:15g}"
            line += f" {tc['D1_1_100'][i]:15g}"
            line += f" {tc0[i]:15g}"
            line += f" {diff[i]:15g}"
            f.write(line + "\n")

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)


def drain_smoothing(xdiff, xrange, newton=False):
    sat = xdiff / xrange
    if not newton:
        f = sat.copy()
    else:
        cof1 = -1.0 / (xrange) ** 3
        cof2 = 2.0 / (xrange) ** 2
        f = cof1 * xdiff**3 + cof2 * xdiff**2
    f[sat < 0] = 0
    f[sat > 1] = 1
    return f


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

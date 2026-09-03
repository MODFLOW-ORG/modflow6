"""
Test MAW SHUT_OFF option.

Two cases share the same Thiem-well setup (rate=-2000, head_limit=-0.4):
  a) HEAD_LIMIT only (baseline) — rate tapers but never reaches exactly 0
  b) HEAD_LIMIT + SHUT_OFF  — well shuts off (rate == 0) when potential flow
     falls below minrate and the well head is below head_limit

The SHUT_OFF keyword takes minrate and maxrate.  minrate is the threshold
below which the well shuts off; maxrate is the threshold above which it
reactivates.  For the shut-off case to trigger, minrate must be chosen large
enough that the taper will cross it during the stress period.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["maw_shutoff_a", "maw_shutoff_b"]

# Grid / solver parameters shared by both cases
nlay, nrow, ncol = 1, 101, 101
nper = 1
perlen = [500.0]
nstp = [50]
tsmult = [1.2]
delr = delc = 142.0
top = 0.0
botm = [-1000.0]
strt = 0.0
hk = 10.0

nouter, ninner = 100, 100
hclose, rclose, relax = 1e-6, 1e-6, 1.0

# MAW parameters
wellbottom = -1000.0
rate = -2000.0
head_limit = -0.4
# SHUT_OFF thresholds: well shuts off once potential flow < minrate
# and reactivates only if it exceeds maxrate.  With the Thiem setup used here,
# the HEAD_LIMIT-throttled rate drops from 2000 to ~616 over 500 days, so
# minrate=900 ensures the shutoff fires before the simulation ends.
shutoff_minrate = 900.0
shutoff_maxrate = 1200.0


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    tdis_rc = [(perlen[0], nstp[0], tsmult[0])]

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws)

    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )

    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(sim.ims, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )

    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=0,
        ss=1.0e-5,
        sy=0.1,
        steady_state={0: False},
        transient={0: True},
        filename=f"{name}.sto",
    )

    # Period data: case a has HEAD_LIMIT only; case b adds SHUT_OFF
    perioddata = {
        0: [
            [0, "rate", rate],
            [0, "head_limit", head_limit],
        ]
    }
    if idx == 1:
        perioddata[0].append([0, "shut_off", shutoff_minrate, shutoff_maxrate])

    mawo_dict = {
        f"{name}.maw.obs.csv": [
            ("m1head", "head", (0,)),
            ("m1rate", "rate", (0,)),
        ]
    }

    flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations=mawo_dict,
        packagedata=[[0, 0.15, wellbottom, strt, "THIEM", 1]],
        connectiondata=[[0, 0, (0, 50, 50), 0.0, wellbottom, 0.0, 0.0]],
        perioddata=perioddata,
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
        printrecord=[("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    fpth = os.path.join(test.workspace, f"{name}.maw.obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except Exception:
        assert False, f'could not load data from "{fpth}"'

    rates = tc["M1RATE"]
    heads = tc["M1HEAD"]

    if idx == 0:
        # HEAD_LIMIT only: rate tapers once the well head reaches head_limit, but
        # never reaches exactly 0 (qpot stays positive as long as h_aq > head_limit).
        assert rates.min() < -shutoff_minrate, (
            "Extraction rate should drop below shutoff_minrate "
            f"({-shutoff_minrate}) with HEAD_LIMIT throttling active; "
            f"min rate = {rates.min()}"
        )
        assert not np.any(np.isclose(rates, 0.0)), (
            "Without SHUT_OFF, the well rate should not snap to exactly 0; "
            f"rates = {rates}"
        )
    else:
        # SHUT_OFF: well must shut off (rate == 0) at some timestep
        assert np.any(np.isclose(rates, 0.0)), (
            "With SHUT_OFF active, the well rate should reach 0 at some timestep "
            f"but min rate = {rates.min()}"
        )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

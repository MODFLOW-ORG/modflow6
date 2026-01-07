# Test NPF highest_cell_saturation option


import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "gwf_npf_dcsa",
    "gwf_npf_dcsb",
    "gwf_npf_dcsc",
]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 1, 4
    nper = 1
    perlen = [1.0]
    nstp = [1]
    tsmult = [1.0]
    delr = 1.0
    delc = 1.0
    top = np.array([[1.0, 1.0, 0.5, 0.5]], dtype=float)
    botm = np.array([[[0.0, 0.0, -1.0, -1.0]]], dtype=float)
    if idx == 0:
        strt = np.array([[[0.25, 0.25, 0.0, 0.0]]], dtype=float)
    else:
        strt = np.array([[[0.25, 0.25, 0.5, 0.5]]], dtype=float)
    hk = 1.0

    c = {0: [((0, 0, 0), strt[0, 0, 0]), ((0, 0, ncol - 1), strt[0, 0, ncol - 1])]}

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = "flow"
    ws = test.workspace

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name="mf6",
        sim_ws=ws,
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
        rcloserecord=500.0,
    )

    if idx in (0, 1):
        newtonoptions = "NEWTON UNDER_RELAXATION"
    else:
        newtonoptions = None

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        newtonoptions=newtonoptions,
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

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        alternative_cell_averaging="AMT-HMK",
        highest_cell_saturation=True,
        icelltype=1,
        k=hk,
        k33=hk,
    )

    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "LAST")],
        printrecord=[("BUDGET", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    ws = test.workspace
    name = "flow"

    if idx in (0, 1):
        sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
        gwf = sim.get_model()

        hds = gwf.output.head().get_data().flatten()
        if idx == 0:
            answer = np.array([0.25, 0.15225742, 0.02386603, 0.0], dtype=float)
        else:
            answer = np.array([0.25, 0.38107678, 0.4666998, 0.5], dtype=float)

        print(f"head: {hds}")
        print(f"answer: {answer}")
        print(f"diff: {hds - answer}")
        assert np.allclose(hds, answer), (
            "simulated head do not match with known solution."
        )
    else:
        buff = test.buffs[0]
        assert any(
            "HIGHEST_CELL_SATURATION option cannot be used when NEWTON option in not"
            in l
            for l in buff
        )


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

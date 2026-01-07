"""
Test the dependent_variable_scale option for gwt for a one-dimensional model grid
of square cells.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

dewatered = [False, True, False, True]
perched = [False, False, True, True]
chd1_lst = (-0.5, 0.0)
cases = [f"vcv{n:02d}_02" for n in range(len(dewatered))]


def build_model(idx, ws, vhfb=False):
    nlay, nrow, ncol = 2, 1, 2
    nper = 2
    perlen = [1.0] * nper
    nstp = [1] * nper
    tsmult = [1.0] * nper
    delr = 1.0
    delc = 1.0
    top = 1.0
    botm = [0.0, -1.0]

    chd0 = 0.5

    strt = chd0

    icelltype = 1
    hk = 100.0  # 50000.0 # 100.0
    area = delc * delr
    vk = 0.001
    if vhfb:
        if dewatered[idx]:
            v0 = 1.0 / ((area * vk) / (0.5 * (chd0 - botm[0])))
            # v1 should be 0.0
            v1 = 0.0
        else:
            v0 = 1.0 / ((area * vk) / (0.5 * (chd0 - botm[0])))
            v1 = 1.0 / ((area * vk) / (0.5 * (botm[0] - botm[1])))
        vcont = 1.0 / (v0 + v1)
    else:
        pass

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="simple",
        outer_dvclose=1e-8,
        outer_maximum=200,
        inner_dvclose=1e-9,
        inner_maximum=100,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        print_input=True,
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="feet",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    if dewatered[idx]:
        cvoptions = "variablecv dewatered"
    else:
        cvoptions = "variablecv"

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        print_flows=True,
        save_flows=True,
        cvoptions=cvoptions,
        perched=perched[idx],
        icelltype=icelltype,
        k=hk,
        k33=vk,
    )

    chd_spd = {}
    for n in range(nper):
        spd = [((0, 0, 0), chd0)]
        spd += [((1, 0, 0), chd1_lst[n])]
        chd_spd[n] = spd
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(spd),
        stress_period_data=chd_spd,
    )

    obs_data = {
        "flowja.obs.csv": [
            ("h0", "head", (0, 0, 1)),
            ("h1", "head", (1, 0, 1)),
            ("vflow", "FLOW-JA-FACE", (0, 0, 1), (1, 0, 1)),
        ]
    }
    obs = flopy.mf6.ModflowUtlobs(gwf, continuous=obs_data)

    return sim


def build_models(idx, test):
    ws = test.workspace
    sim = build_model(idx, ws, vhfb=True)

    return sim, None


def check_results(idx, test):
    ws = test.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws, verbosity_level=0)
    gwf = sim.get_model()

    nper = sim.tdis.nper.array
    delr = gwf.dis.delr.array[1]
    delc = gwf.dis.delc.array[0]
    area = delr * delc

    cellids = [(0, 0, 1), (1, 0, 1)]
    botm0 = gwf.dis.botm.array[cellids[0]]
    botm1 = gwf.dis.botm.array[cellids[1]]
    vk = gwf.npf.k33.array[cellids[0]]

    obs = gwf.obs.output.obs().get_data()
    h0 = obs["H0"]
    h1 = obs["H1"]

    answer = obs["VFLOW"]

    is_dewatered = dewatered[idx]
    is_perched = perched[idx]

    vflow = []
    for n in range(nper):
        if is_dewatered and h1[n] < botm0:
            v0 = 1.0 / ((area * vk) / (0.5 * (h0[n] - botm0)))
            v1 = 0.0
        else:
            v0 = 1.0 / ((area * vk) / (0.5 * (h0[n] - botm0)))
            v1 = 1.0 / ((area * vk) / (0.5 * (botm0 - botm1)))
        vcont = 1.0 / (v0 + v1)
        vflow.append(vcont * (h1[n] - h0[n]))

        if is_perched and h1[n] < botm0:
            qcorr = vcont * (h1[n] - botm0)
            vflow[n] -= qcorr

    vflow = np.array(vflow)

    print(f"H0: {h0}")
    print(f"H1: {h1}")
    print(f"Vertical flow: {answer}")
    print(f"Answer:        {vflow}")

    assert np.allclose(vflow, answer), "simulated results not equal to the answer"


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_results(idx, t),
    )
    test.run()

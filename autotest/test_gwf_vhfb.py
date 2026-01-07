"""
Test the dependent_variable_scale option for gwt for a one-dimensional model grid
of square cells.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cell_averaging = [None, None, None, None] + ["logarithmic"] * 2
xt3d = [False, False, True, True] + [False] * 2
newton = [False, True, False, True] + [False, True]
cases = [f"vhfb{n:02d}" for n in range(len(xt3d))]


def build_model(idx, ws, vhfb=False):
    nlay, nrow, ncol = 2, 11, 11
    nper = 1
    perlen = [2.0]
    nstp = [1]
    tsmult = [1.0]
    delr = 100.0
    delc = 100.0
    top = 100.0
    botm = [50.0, 0.0]
    strt = 75.0
    hk = 100.0  # 50000.0
    if vhfb:
        vk = None
        khfb = 0.001
    else:
        vk = 0.001
    width = 50.0
    icelltype = 1

    if newton[idx]:
        linear_acceleration = "bicgstab"
        newtonoptions = "newton"
    else:
        linear_acceleration = "cg"
        newtonoptions = None

    if xt3d[idx] and not newton[idx]:
        linear_acceleration = "bicgstab"

    w = {
        0: [
            ((1, 5, 5), -1000.0),
        ]
    }

    chd0 = 100.0
    chd1 = 75.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        print_input=True,
        newtonoptions=newtonoptions,
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="simple",
        linear_acceleration=linear_acceleration,
        outer_dvclose=1e-8,
        outer_maximum=200,
        inner_dvclose=1e-9,
        inner_maximum=100,
    )
    sim.register_ims_package(ims, [gwf.name])

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

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        alternative_cell_averaging=cell_averaging[idx],
        xt3doptions=xt3d[idx],
        icelltype=icelltype,
        k=hk,
        k33=vk,
    )

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(w),
        stress_period_data=w,
        save_flows=False,
        pname="WEL-1",
    )

    spd = [((0, i, 0), chd0) for i in range(nrow)]
    spd += [((0, i, ncol - 1), chd1) for i in range(nrow)]
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        maxbound=len(spd),
        stress_period_data=spd,
    )

    if vhfb:
        spd = []
        for i in range(nrow):
            spd += [((0, i, j), (1, i, j), khfb / width) for j in range(ncol)]

        hfb = flopy.mf6.ModflowGwfhfb(
            gwf,
            maxhfb=len(spd),
            stress_period_data=spd,
        )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "ALL"), ("HEAD", "ALL")],
    )

    return sim


def build_models(idx, test):
    ws = test.workspace
    sim = build_model(idx, ws, vhfb=True)

    ws = test.workspace / "mf6"
    mc = build_model(idx, ws, vhfb=False)

    return sim, mc


def check_results(test):
    ws = test.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwf = sim.get_model()
    flow = gwf.output.budget().get_data(text="FLOW-JA-FACE")[0]

    ws = test.workspace / "mf6"
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwf = sim.get_model()
    answer = gwf.output.budget().get_data(text="FLOW-JA-FACE")[0]
    diff = flow - answer
    diffmax = np.abs(diff).max()

    print(f"flow:       {flow.flatten()}\n")
    print(f"answer:     {answer.flatten()}\n")
    print(f"difference: {diff.flatten()}\n")

    assert diffmax < 0.0001, "FLOW-JA-FACE flows are not close to the defined answer."


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        compare="mf6",
        build=lambda t: build_models(idx, t),
        check=lambda t: check_results(t),
    )
    test.run()

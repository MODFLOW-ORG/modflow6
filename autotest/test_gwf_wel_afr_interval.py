"""
AUTO_FLOW_REDUCE must not reduce a withdrawal from a confined cell.

The NPF Package now carries BOTREDUCE and TOPREDUCE, the elevation bounds of
the interval over which AUTO_FLOW_REDUCE tapers a withdrawal. They are empty
unless a flow formulation registers and narrows them, and WEL falls back to
the cell bottom with the convertible-cell rule when they are empty. That
fallback changed the test WEL applies, from

    if (icelltype == 0) cycle

to

    if (icelltype == 0 .and. .not. associated(botreduce)) cycle

so a confined cell is still skipped when no formulation is active. The
existing WEL tests all use convertible cells, so this covers the confined
case the condition now depends on.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["welafr-confined"]

top, botm = 1.0, -40.0
hchd = 1.0
qrequest = -0.1
k = 0.01
afr = 0.1


def build_models(idx, test):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(400.0, 40, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=200,
        inner_maximum=200,
        outer_dvclose=1.0e-9,
        inner_dvclose=1.0e-10,
        linear_acceleration="bicgstab",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=2, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=hchd)
    # confined: the reduction must not apply, with or without an interval
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, sy=0.2, transient=True)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), hchd]])
    flopy.mf6.ModflowGwfwel(
        gwf, auto_flow_reduce=afr, stress_period_data=[[(0, 0, 1), qrequest]]
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budget_filerecord=name + ".bud",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, None


def check_output(idx, test):
    gwf = test.sims[0].gwf[0]
    cbc = gwf.output.budget()
    hobj = gwf.output.head()

    worst = 0.0
    for t in cbc.get_times():
        rec = cbc.get_data(text="WEL", totim=t)[0]
        q = float(np.asarray(rec["q"]).ravel()[0])
        worst = max(worst, abs(q - qrequest))
        assert abs(q - qrequest) < 1.0e-9, (
            f"a confined cell must deliver the requested rate {qrequest}, "
            f"got {q} at time {t}"
        )

    hlast = hobj.get_data(totim=cbc.get_times()[-1]).flatten()[1]
    print(
        f"confined cell: rate held at {qrequest} (max deviation {worst:.1e}); "
        f"final head {hlast:.4f}"
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

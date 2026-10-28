"""
AUTO_FLOW_REDUCE with a flow-reduction provider attached, or not.

WEL asks an attached provider (a flow formulation such as SWI) how much fluid
a cell holds before it tapers a withdrawal. With no provider the taper is the
one it always had: measured up from the cell bottom, over the fraction or
length the user gave, and applied to convertible cells only. Two things that
rule has to keep doing:

welafr-confined
    A confined cell is not reduced when no provider is attached. The
    convertible-cell test used to be `if (icelltype == 0) cycle`; it now
    depends on whether a provider narrowed the interval, so the plain confined
    case is covered here.

welafr-longtaper
    A global FLOW_REDUCTION_LENGTH larger than the cell thickness puts the top
    of the taper above the cell top. The taper must still run to bottom plus
    length: the well is partly reduced while the head is inside the taper and
    returns to the full rate once the head is above it. An early version of
    the provider code capped the head at the cell top even with no provider,
    which froze the rate at the taper value for the cell top (0.5 here)
    however high the head rose. The case is the one reported in the review of
    the pull request: a 10 m cell, a 20 m length, heads of 15 and 25 m.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["welafr-confined", "welafr-longtaper"]

# welafr-confined
top, botm = 1.0, -40.0
hchd = 1.0
qrequest = -0.1
k = 0.01
afr = 0.1

# welafr-longtaper
lt_top, lt_botm = 10.0, 0.0
lt_heads = [15.0, 25.0]
lt_length = 20.0
lt_q = -1.0
lt_k = 10.0


def taper(tp, bt, x):
    """sQSaturation of the SmoothingModule: the cubic taper of AUTO_FLOW_REDUCE."""
    s = (x - bt) / (tp - bt)
    if s < 0.0:
        return 0.0
    if s < 1.0:
        return -2.0 * s**3 + 3.0 * s**2
    return 1.0


def build_confined(sim, name):
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
    return gwf


def build_longtaper(sim, name):
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=1.0, delc=1.0, top=lt_top, botm=lt_botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=lt_heads[0])
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=lt_k)
    # both ends held at the period's head; the well in the middle cell
    chd = {kper: [[(0, 0, 0), h], [(0, 0, 2), h]] for kper, h in enumerate(lt_heads)}
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfwel(
        gwf,
        auto_flow_reduce=lt_length,
        flow_reduction_length=True,
        stress_period_data=[[(0, 0, 1), lt_q]],
    )
    return gwf


def build_models(idx, test):
    name = "mymodel"
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    if cases[idx] == "welafr-confined":
        flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(400.0, 40, 1.0)])
    else:
        flopy.mf6.ModflowTdis(
            sim, nper=len(lt_heads), perioddata=len(lt_heads) * [(1.0, 1, 1.0)]
        )
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_maximum=200,
        inner_maximum=200,
        outer_dvclose=1.0e-9,
        inner_dvclose=1.0e-10,
        linear_acceleration="bicgstab",
    )
    if cases[idx] == "welafr-confined":
        gwf = build_confined(sim, name)
    else:
        gwf = build_longtaper(sim, name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=name + ".hds",
        budget_filerecord=name + ".bud",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, None


def check_confined(test):
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


def check_longtaper(test):
    gwf = test.sims[0].gwf[0]
    cbc = gwf.output.budget()
    hobj = gwf.output.head()
    for t, hend in zip(cbc.get_times(), lt_heads):
        rec = cbc.get_data(text="WEL", totim=t)[0]
        q = float(np.asarray(rec["q"]).ravel()[0])
        h = hobj.get_data(totim=t).flatten()[1]
        expected = lt_q * taper(lt_botm + lt_length, lt_botm, h)
        print(
            f"ends at {hend}: well-cell head {h:.4f}, rate {q:.6f}, "
            f"taper to {lt_botm + lt_length} gives {expected:.6f}"
        )
        assert abs(q - expected) < 1.0e-8, (
            f"the taper must run to bottom plus length; got {q}, expected {expected}"
        )
    # above the taper the well pumps at the requested rate
    assert abs(q - lt_q) < 1.0e-9, f"full rate expected above the taper, got {q}"


def check_output(idx, test):
    if cases[idx] == "welafr-confined":
        check_confined(test)
    else:
        check_longtaper(test)


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

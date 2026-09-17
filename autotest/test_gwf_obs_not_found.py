"""
Test the warning issued when nothing is observed at an observation location.

An observation on a traditional stress package identifies its location either
by a boundname or by a cellid. Either one can fail to match any boundary in
the package: the boundname may be misspelled or unused, and the cell may
simply not contain a boundary. The observation is not an error, because the
boundary list of a stress package is allowed to change between stress periods,
but nothing is observed at that location and the DNODATA value is reported
instead of a rate. MODFLOW 6 issues a warning so that the empty observation is
not mistaken for a simulated value of zero.

The model is a one-layer strip with constant heads at both ends and a single
well, and it carries four WEL Package observations:

* found_by_name   - the boundname of the well, which must not be warned about
* missing_by_name - a boundname no well uses, which must be warned about
* found_by_cell   - the cell containing the well, which must not be warned
                    about
* missing_by_cell - a cell with no boundary in it, which must be warned about

The warning names the observation and the package, and reports the boundname
or the cell that was not matched, so the warnings are matched against the
observations that are expected to produce them.

Two stress periods are simulated. Nothing about the well changes between them,
so an unmatched observation is unmatched in both, and the warning must still
be reported only once rather than once per stress period.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["obsnotfound"]

# grid and aquifer properties
nlay, nrow, ncol = 1, 1, 6
delr = delc = 10.0
top = 10.0
botm = 0.0
hk = 10.0
strt = top

# the single well, its boundname and its cell
wellcol = 2
wellname = "well_a"
wellrate = -1.0

# a cell with no boundary in it; not one of the constant-head cells
emptycol = 3

# solver settings
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-6, 0.97

# observations: (name, id string passed to the observation, expect a warning,
# the kind of location the warning must report)
observations = [
    ("found_by_name", wellname, False, None),
    ("missing_by_name", "no_such_well", True, "boundname"),
    ("found_by_cell", (0, 0, wellcol), False, None),
    ("missing_by_cell", (0, 0, emptycol), True, "cell"),
]

# DNODATA, the value reported when nothing is observed
DNODATA = 3.0e30


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", sim_ws=test.workspace)
    # two stress periods, so an unmatched observation is searched for twice
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=2, perioddata=[(1.0, 1, 1.0), (1.0, 1, 1.0)]
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        relaxation_factor=relax,
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, save_flows=True)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), strt], [(0, 0, ncol - 1), strt]],
    )

    # the well is present in both stress periods, so the observations that do
    # match it must not be warned about for either period
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        pname="wel",
        boundnames=True,
        print_input=True,
        print_flows=True,
        save_flows=True,
        stress_period_data={
            0: [[(0, 0, wellcol), wellrate, wellname]],
            1: [[(0, 0, wellcol), wellrate, wellname]],
        },
    )
    wel.obs.initialize(
        filename=f"{name}.wel.obs",
        continuous={
            f"{name}.wel.obs.csv": [
                (obsname, "wel", obsid) for obsname, obsid, _, _ in observations
            ]
        },
    )
    return sim


def reported_warnings(ws):
    """Warnings about unmatched observations, from mfsim.lst.

    Returns a list of (kind, location, observation name, package name), where
    kind is "boundname" or "cell".
    """
    with open(os.path.join(ws, "mfsim.lst"), "r") as f:
        # the warnings are wrapped, so the line breaks are removed before the
        # message is searched for
        text = " ".join(f.read().split())
    pattern = (
        r'No boundary matching (boundname|cell) "([^"]*)" was found '
        r'for observation "([^"]*)" in package "([^"]*)"\. '
        r"The DNODATA value will be returned for this observation\."
    )
    return [m.groups() for m in re.finditer(pattern, text)]


def check_output(idx, test):
    name = cases[idx]
    warnings = reported_warnings(test.workspace)

    # every observation that matches nothing must be warned about exactly
    # once, and no other observation may be warned about; counting rather
    # than collecting into a set also checks that the warning is deduplicated
    # across the two stress periods
    warned = [obsname.lower() for _, _, obsname, _ in warnings]
    expected = sorted(obsname for obsname, _, expect, _ in observations if expect)
    assert sorted(warned) == expected, (
        f"warnings were reported for {sorted(warned)}, expected {expected}"
    )

    # the warning must report the kind of location that was not matched and
    # name the package the observation belongs to
    kinds = {obsname.lower(): kind for kind, _, obsname, _ in warnings}
    packages = {obsname.lower(): pkg.lower() for _, _, obsname, pkg in warnings}
    for obsname, _, expect, kind in observations:
        if not expect:
            continue
        assert kinds[obsname] == kind, (
            f"observation {obsname} was warned about as a {kinds[obsname]}, "
            f"expected {kind}"
        )
        assert packages[obsname] == "wel", (
            f"observation {obsname} was attributed to package "
            f"{packages[obsname]}, expected wel"
        )

    # the location that was not matched must be the one that was requested,
    # and must be reported without the padding the id string carries
    located = {obsname.lower(): loc.lower() for _, loc, obsname, _ in warnings}
    assert located["missing_by_name"] == "no_such_well", (
        f"the warning reported boundname {located['missing_by_name']}, "
        "expected no_such_well"
    )
    cellid = " ".join(str(i + 1) for i in (0, 0, emptycol))
    assert located["missing_by_cell"] == cellid, (
        f"the warning reported cell {located['missing_by_cell']}, expected {cellid}"
    )

    # the warning claims DNODATA is returned, so check that it is, and that
    # the observations that do match the well report its rate instead
    obs = flopy.utils.Mf6Obs(
        os.path.join(test.workspace, f"{name}.wel.obs.csv"), isBinary=False
    ).get_data()
    for obsname, _, expect, _ in observations:
        values = obs[obsname.upper()]
        if expect:
            assert np.all(values == DNODATA), (
                f"observation {obsname} reported {values}, expected DNODATA"
            )
        else:
            assert np.allclose(values, wellrate), (
                f"observation {obsname} reported {values}, expected {wellrate}"
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

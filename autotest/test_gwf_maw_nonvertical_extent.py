"""
Test the warning issued when a non-vertical multi-aquifer well connection spans
a horizontal distance greater than the cell it is connected to.

A connection listed in the ANGLEDATA block spans a horizontal distance of
LW * sin(ANGLE), where LW is the in-cell screen length. The screen cannot be
longer than the connected cell, so a connection that spans more than the cell
is given a saturated conductance that is too large. MODFLOW 6 compares the
horizontal distance to the maximum horizontal extent of the connected cell and
issues a warning when it is exceeded.

The cell extent is calculated from the cell vertices where they are defined and
is estimated from the cell area where they are not. The vertices are optional
for a DISU grid, so it is run both with and without them. All four grids in this
test are the same square cells, so all four must report the same extent, and
only the DISU grid without vertices must report it as estimated.

The same four wells are simulated on each grid:

* inside      - a 60 degree connection with a screen that spans the cell
                thickness; the horizontal distance is less than the cell extent
                and must not be warned about
* outside     - an 85 degree connection with a screen that spans the cell
                thickness; the horizontal distance is greater than the cell
                extent and must be warned about
* lat_inside  - a horizontal connection with a CONN_LENGTH less than the cell
                extent, which must not be warned about
* lat_outside - a horizontal connection with a CONN_LENGTH greater than the cell
                extent, which must be warned about
* specified   - the same as lat_outside but using the SPECIFIED conductance
                equation, for which the program does not calculate the
                conductance and must not report it as too large
"""

import math
import os
import re

import flopy
import numpy as np
import pytest
from flopy.utils.gridutil import get_disu_kwargs, get_disv_kwargs
from framework import TestFramework

cases = ["mawnvext_dis", "mawnvext_disv", "mawnvext_disu", "mawnvext_disuv"]
grids = ["dis", "disv", "disu", "disu_vertices"]

# grid and aquifer properties; square cells so the maximum horizontal extent of
# a cell is the diagonal, and is the same for all three grids
nlay, nrow, ncol = 1, 1, 11
delr = delc = 20.0
top = 10.0
botm = 0.0
hk = 10.0
strt = top
extent = math.sqrt(2.0) * delr

# multi-aquifer well properties
radius = 0.25
sradius = 0.5
hks = hk
satcond = 1.0
mawrate = -1.0

# wells: (column, angle in degrees, connection length, conductance equation,
# expect a warning)
wells = [
    (1, 60.0, None, "MEAN", False),
    (3, 85.0, None, "MEAN", True),
    (5, 90.0, 15.0, "MEAN", False),
    (7, 90.0, 50.0, "MEAN", True),
    (9, 90.0, 50.0, "SPECIFIED", True),
]

# the conductance clause of the warning depends on the conductance equation
CALCULATED = "The calculated saturated conductance is correspondingly too large."
SPECIFIED = (
    "The specified saturated conductance is applied to a screen "
    "that is longer than the cell."
)

# a horizontal connection is saturated over the well diameter
lat_top = 0.5 * (top + botm) + radius
lat_bot = 0.5 * (top + botm) - radius

nouter, ninner = 100, 100
hclose, rclose = 1e-9, 1e-9


def horizontal_distance(angle, conn_len):
    """Horizontal distance spanned by a connection."""
    omega = math.radians(angle)
    if conn_len is None:
        lw = (top - botm - 2.0 * radius * math.sin(omega)) / math.cos(omega)
    else:
        lw = conn_len
    return lw * math.sin(omega)


def add_grid(gwf, grid):
    """Add the discretization package, and return a cellid for each column."""
    if grid == "dis":
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
        return [(0, 0, j) for j in range(ncol)]
    elif grid == "disv":
        kwargs = get_disv_kwargs(nlay, nrow, ncol, delr, delc, top, [botm])
        flopy.mf6.ModflowGwfdisv(gwf, **kwargs)
        return [(0, j) for j in range(ncol)]
    else:
        # the cell vertices are optional; the extent is estimated from the cell
        # area when they are not given
        kwargs = get_disu_kwargs(
            nlay,
            nrow,
            ncol,
            delr,
            delc,
            top,
            [botm],
            return_vertices=grid == "disu_vertices",
        )
        flopy.mf6.ModflowGwfdisu(gwf, **kwargs)
        return [j for j in range(ncol)]


def build_models(idx, test):
    name = cases[idx]
    grid = grids[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    cellids = add_grid(gwf, grid)

    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    # confined aquifer so the screen saturation is 1.0 and the full saturated
    # conductance is applied
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk, save_flows=True)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[cellids[0], strt], [cellids[-1], strt]],
    )

    packagedata = []
    connectiondata = []
    angledata = []
    perioddata = []
    for ifno, (col, angle, conn_len, condeqn, _) in enumerate(wells):
        packagedata.append([ifno, radius, botm, strt, condeqn, 1])
        if conn_len is None:
            scrn_top, scrn_bot = top, botm
        else:
            scrn_top, scrn_bot = lat_top, lat_bot
        connectiondata.append([ifno, 0, cellids[col], scrn_top, scrn_bot, hks, sradius])
        if conn_len is None:
            angledata.append([ifno, 0, angle])
        else:
            angledata.append([ifno, 0, angle, conn_len])
        perioddata.append([ifno, "rate", mawrate])

    flopy.mf6.ModflowGwfmaw(
        gwf,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        non_vertical_wells=True,
        packagedata=packagedata,
        connectiondata=connectiondata,
        angledata=angledata,
        perioddata={0: perioddata},
    )
    return sim


def warned_connections(ws):
    """Wells warned about, and the extent reported for each, from mfsim.lst."""
    with open(os.path.join(ws, "mfsim.lst"), "r") as f:
        # the warnings are wrapped, so the line breaks are removed before the
        # message is searched for
        text = " ".join(f.read().split())
    pattern = (
        r"The horizontal distance spanned by maw well (\d+) connection \d+ "
        r"\(([-+.0-9eE]+)\) is greater than the (estimated )?maximum "
        r"horizontal extent of the connected cell \(([-+.0-9eE]+)\), so the "
        r"screen extends beyond the cell it is connected to\. (.+?) Reduce "
        r"ANGLE"
    )
    return {
        int(m.group(1)): (
            float(m.group(2)),
            float(m.group(4)),
            m.group(5),
            m.group(3) is not None,
        )
        for m in re.finditer(pattern, text)
    }


def check_output(idx, test):
    warned = warned_connections(test.workspace)
    # only the disu grid without vertices estimates the extent from the area
    estimate = grids[idx] == "disu"
    expected = {i + 1 for i, w in enumerate(wells) if w[4]}
    assert set(warned) == expected, (
        f"warned about wells {sorted(warned)}, expected {sorted(expected)}"
    )

    for ifno, (col, angle, conn_len, condeqn, warn) in enumerate(wells):
        if not warn:
            continue
        hlen, reported, clause, estimated = warned[ifno + 1]
        assert np.isclose(hlen, horizontal_distance(angle, conn_len), rtol=1e-6), (
            f"well {ifno + 1} horizontal distance {hlen} is not the expected "
            f"{horizontal_distance(angle, conn_len)}"
        )
        # the extent is calculated from the cell vertices where they are
        # defined and estimated from the cell area where they are not, and
        # must be the same for all four grids
        assert np.isclose(reported, extent, rtol=1e-6), (
            f"well {ifno + 1} cell extent {reported} is not the expected {extent}"
        )
        assert estimated == estimate, (
            f"well {ifno + 1} reported the extent as "
            f"{'estimated' if estimated else 'calculated'} on the "
            f"{grids[idx]} grid"
        )
        # the program does not calculate the conductance of a SPECIFIED
        # connection, so it must not report it as too large
        expect = SPECIFIED if condeqn == "SPECIFIED" else CALCULATED
        assert clause == expect, (
            f"well {ifno + 1} reported '{clause}', expected '{expect}'"
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

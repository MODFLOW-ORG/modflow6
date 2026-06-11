"""
Test that TVK time-series K values are correctly applied in stress periods
where the TVK package has no explicit period data (reuse periods).

When a TVK file has a PERIOD block only for period 1, the K values for
periods 2+ come from time-series interpolation.  This test verifies that
TvBase.ad re-applies those IDM-updated values at every time step, including
time steps that fall inside reuse periods.

Model layout
------------
3 layers, 1 row, 1 column.  CHD at the top cell (h=1.0) and bottom cell
(h=0.0).  K11=K33=1.0 everywhere initially.  TVK varies K33 of the bottom
cell only, via a linear time series K33(t) = 1 + t.

Three stress periods of 1 day each, 1 step each.
TVK has a PERIOD 1 block only; periods 2-3 are reuse periods.

Analytical head at the middle cell for serial 1-D vertical flow:

    h_mid = (K + 1) / (3*K + 1)

where K = K33 of the bottom cell.  MODFLOW 6 evaluates time series at the
midpoint of each time step, so for 1-day steps K33(t_mid) = 1 + t_mid:

    period 1  (t_mid=0.5, K33=1.5)  h_mid = 5/11
    period 2  (t_mid=1.5, K33=2.5)  h_mid = 7/17
    period 3  (t_mid=2.5, K33=3.5)  h_mid = 9/23
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["tvk06ts"]

nlay, nrow, ncol = 3, 1, 1
nper = 3
perlen = [1.0] * nper
nstp = [1] * nper
total_time = float(sum(perlen))

cellid_top = (0, 0, 0)
cellid_bot = (2, 0, 0)
cellid_mid = (1, 0, 0)


def _expected_head(k33):
    """Analytical middle-cell head for a uniform-K33 serial column."""
    return (k33 + 1.0) / (3.0 * k33 + 1.0)


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]
    gwfname = f"gwf_{name}"

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=[(p, s, 1.0) for p, s in zip(perlen, nstp)],
    )

    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=gwfname,
        model_nam_file=f"{gwfname}.nam",
    )
    imsgwf = flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", filename=f"{gwfname}.ims"
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=0.0,
        botm=[-1.0, -2.0, -3.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.5)

    tvk_filename = f"{gwfname}.npf.tvk"
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0, k33=1.0)

    # TVK: PERIOD 1 only; K33 of bottom cell via linear time series K33(t)=1+t
    tvkspd = {0: [[cellid_bot, "K33", "k33_series"]]}
    tvk = flopy.mf6.ModflowUtltvk(
        npf,
        print_input=True,
        perioddata=tvkspd,
        filename=tvk_filename,
    )
    tvk.ts.initialize(
        filename=f"{gwfname}.npf.tvk.ts",
        timeseries=[(float(t), 1.0 + float(t)) for t in range(int(total_time) + 1)],
        time_series_namerecord="k33_series",
        interpolation_methodrecord="linear",
    )

    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[cellid_top, 1.0], [cellid_bot, 0.0]],
        pname="CHD-1",
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST")],
    )

    return sim, None


def check_output(idx, test):
    """
    Verify that the middle-cell head matches the analytical value at each
    stress period, confirming that the TVK time-series K33 value is applied
    in reuse stress periods.
    """
    gwfname = f"gwf_{cases[idx]}"
    fpth = os.path.join(test.workspace, f"{gwfname}.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    head = hobj.get_alldata()

    assert len(head) == nper, f"Expected {nper} head records, got {len(head)}"

    # K33 is evaluated at the midpoint of each 1-day time step: t_mid = k + 0.5
    k33_per_period = [1.0 + float(k) + 0.5 for k in range(nper)]
    expected = [_expected_head(k33) for k33 in k33_per_period]

    for kper, (h_exp, k33) in enumerate(zip(expected, k33_per_period)):
        h = float(head[kper, cellid_mid[0], cellid_mid[1], cellid_mid[2]])
        assert np.isclose(h, h_exp, atol=1e-6), (
            f"Period {kper + 1}: middle-cell head {h:.8f} does not match "
            f"analytical value {h_exp:.8f} for K33={k33:.1f}.  "
            f"TVK time-series values are not being applied in reuse stress "
            f"periods."
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

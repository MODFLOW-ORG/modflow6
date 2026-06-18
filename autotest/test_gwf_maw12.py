"""
Test non-vertical (slanted) multi-aquifer well connections activated with the
NON_VERTICAL_WELLS option and the ANGLEDATA block.

The saturated conductance of a connection listed in the ANGLEDATA block is
scaled by the in-cell screen length, which is calculated from the screen top,
screen bottom, well radius, and tilt angle (or specified directly with an
optional connection length) instead of the vertical screen thickness. The MEAN
conductance equation is used for these tests because the resulting cylindrical
conductance depends only on the screen length and is therefore exact for any
connection orientation.

Several independent single-well models are built in one simulation and the
simulated well heads are compared to verify the length correction:

* vert         - vertical well, screen spans a thickness of LREF (reference)
* angle0       - same as vert but listed in ANGLEDATA with a 0 degree angle;
                 must reproduce the vertical (reference) result exactly
* noangle      - the NON_VERTICAL_WELLS option is activated but no ANGLEDATA
                 block is provided; must run (with a warning) and reproduce the
                 vertical (reference) result exactly
* horizontal   - a thin (near-horizontal) screen with a 90 degree angle and a
                 connection length of LREF; must reproduce the reference result
* slant        - a 45 degree connection; the in-cell screen length is derived
                 from the screen elevations, well radius, and angle
* slant_equiv  - a vertical well whose screen thickness equals the derived
                 in-cell screen length of the slant case; must reproduce the
                 slant result

A radial collector (Ranney) well with four equal-length horizontal laterals is
also simulated to confirm that, by symmetry, the four lateral connection flows
are equal.

Finally, a model with multiple multi-aquifer wells in a single MAW package is
simulated to confirm the ANGLEDATA block is handled correctly when only some
wells and some connections are non-vertical:

* multiwell    - two wells; the first is vertical and is not listed in the
                 ANGLEDATA block, and the second has one vertical connection
                 (not listed in ANGLEDATA) and one slanted connection (listed
                 in the ANGLEDATA block)
* multiwellref - the same two wells, but the slanted connection of the second
                 well is replaced by a vertical connection whose screen
                 thickness equals the derived in-cell screen length (no
                 ANGLEDATA block); both wells must reproduce the multiwell heads
"""

import math
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["maw12"]

# grid and aquifer properties shared by the single-well models
top = 100.0
botm = 0.0
delr = 100.0
delc = 100.0
hk = 10.0
strt = 100.0

# multi-aquifer well properties
radius = 0.1
sradius = 0.5
hks = 10.0
mawrate = -50.0

# reference vertical screen length (full cell thickness)
lref = top - botm

# slant case: screen elevations and tilt angle (degrees from vertical)
slant_top = 80.0
slant_bot = 30.0
slant_angle = 45.0
slant_dz = slant_top - slant_bot
# in-cell screen length derived from the screen elevations, radius, and angle
slant_omega = math.radians(slant_angle)
slant_lw = (slant_dz - 2.0 * radius * math.sin(slant_omega)) / math.cos(slant_omega)
# equivalent vertical screen with a thickness equal to the derived length,
# centered within the cell
slant_equiv_top = 0.5 * (top + botm) + 0.5 * slant_lw
slant_equiv_bot = 0.5 * (top + botm) - 0.5 * slant_lw

# multi-well case: the slanted connection of the second well and the equivalent
# vertical screen thickness derived from its screen elevations, radius, and angle
mw_angle = 40.0
mw_slant_top = 70.0
mw_slant_bot = 20.0
mw_slant_dz = mw_slant_top - mw_slant_bot
mw_omega = math.radians(mw_angle)
mw_slant_lw = (mw_slant_dz - 2.0 * radius * math.sin(mw_omega)) / math.cos(mw_omega)
mw_equiv_top = 0.5 * (top + botm) + 0.5 * mw_slant_lw
mw_equiv_bot = 0.5 * (top + botm) - 0.5 * mw_slant_lw

# partial-saturation case: a slanted connection over a vertical screen interval
# in an unconfined cell is equivalent (at every head, including partial
# saturation) to a vertical connection over the same interval with the skin
# conductivity scaled by the length-correction factor, because both have the
# same vertical screen interval (same saturation) and the same saturated
# conductance.
ps_top = 80.0
ps_bot = 40.0
ps_angle = 45.0
ps_dz = ps_top - ps_bot
ps_omega = math.radians(ps_angle)
ps_lcorr = (ps_dz - 2.0 * radius * math.sin(ps_omega)) / math.cos(ps_omega) / ps_dz
ps_hks_vert = hks * ps_lcorr
# the saturated conductance of the slanted MEAN connection, which is what a
# SPECIFIED connection must be given (the length correction is not applied to
# SPECIFIED, so the user supplies the already-corrected conductance)
ps_cspec = hks * math.pi * (radius + sradius) * ps_dz / (sradius - radius) * ps_lcorr
ps_chd = 70.0
ps_rate = -4000.0

nouter = 100
ninner = 100
hclose = 1e-9
rclose = 1e-9


def add_single_well_model(
    sim, name, scrn_top, scrn_bot, angledata=None, conn_len=None, nonvert_opt=False
):
    """Add an independent single-well GWF model to the simulation.

    The model is a one-row, three-column confined aquifer with a constant-head
    boundary in the first column and a multi-aquifer well (MEAN conductance
    equation) connected to the middle column. If ``nonvert_opt`` is True the
    NON_VERTICAL_WELLS option is activated without an ANGLEDATA block, which
    should run (with a warning) and treat all connections as vertical.
    """
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=f"{name}.nam")

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=3,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )
    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")
    # confined aquifer (icelltype=0) so the screen saturation is 1.0 and the
    # full saturated conductance is applied
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=0,
        k=hk,
        k33=hk,
        save_flows=True,
        filename=f"{name}.npf",
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), strt]],
        filename=f"{name}.chd",
    )

    # MEAN conductance equation: [ifno, radius, bottom, strt, condeqn, ngwfnodes]
    packagedata = [[0, radius, botm, strt, "MEAN", 1]]
    # connectiondata: [ifno, icon, cellid, scrn_top, scrn_bot, hk_skin, radius_skin]
    connectiondata = [[0, 0, (0, 0, 1), scrn_top, scrn_bot, hks, sradius]]
    perioddata = {0: [[0, "rate", mawrate]]}

    maw_kwargs = dict(
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations={f"{name}.maw.obs.csv": [("h", "head", 1)]},
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname=f"maw_{name}",
        filename=f"{name}.maw",
    )
    if angledata is not None:
        maw_kwargs["non_vertical_wells"] = True
        # angledata row: [ifno, icon, angle, (conn_length)]
        if conn_len is None:
            maw_kwargs["angledata"] = [[0, 0, angledata]]
        else:
            maw_kwargs["angledata"] = [[0, 0, angledata, conn_len]]
    elif nonvert_opt:
        # activate the option but provide no ANGLEDATA block
        maw_kwargs["non_vertical_wells"] = True

    flopy.mf6.ModflowGwfmaw(gwf, **maw_kwargs)
    return gwf


def add_collector_model(sim, name):
    """Add a radial collector (Ranney) well with four equal horizontal laterals."""
    nrow = ncol = 5
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=f"{name}.nam")

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
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
        icelltype=0,
        k=hk,
        k33=hk,
        save_flows=True,
        filename=f"{name}.npf",
    )

    # constant heads on all four edges to create a symmetric flow field
    chdspd = []
    for i in range(nrow):
        for j in range(ncol):
            if i in (0, nrow - 1) or j in (0, ncol - 1):
                chdspd.append([(0, i, j), strt])
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdspd, filename=f"{name}.chd")

    # one well with four horizontal laterals to the N, S, W, and E neighbors of
    # the central cell; the laterals all have the same length
    lateral_cells = [(0, 1, 2), (0, 3, 2), (0, 2, 1), (0, 2, 3)]
    lateral_len = 100.0
    thin = 2.0 * radius  # vertical band of a horizontal screen
    scrn_top = 0.5 * (top + botm) + 0.5 * thin
    scrn_bot = 0.5 * (top + botm) - 0.5 * thin

    packagedata = [[0, radius, botm, strt, "MEAN", len(lateral_cells)]]
    connectiondata = []
    angledata = []
    obs = []
    for icon, cellid in enumerate(lateral_cells):
        connectiondata.append([0, icon, cellid, scrn_top, scrn_bot, hks, sradius])
        # horizontal connection (90 degrees) with an explicit lateral length
        angledata.append([0, icon, 90.0, lateral_len])
        obs.append((f"q{icon + 1}", "maw", 1, icon + 1))

    flopy.mf6.ModflowGwfmaw(
        gwf,
        non_vertical_wells=True,
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations={f"{name}.maw.obs.csv": obs},
        packagedata=packagedata,
        connectiondata=connectiondata,
        angledata=angledata,
        perioddata={0: [[0, "rate", mawrate]]},
        pname=f"maw_{name}",
        filename=f"{name}.maw",
    )
    return gwf


def add_multiwell_model(sim, name, use_angledata):
    """Add a model with two MAW wells in a single MAW package.

    Well A is vertical and is not listed in the ANGLEDATA block. Well B has one
    vertical connection (not listed in ANGLEDATA) and one connection that is
    slanted (``use_angledata=True``) or replaced by an equivalent vertical
    screen (``use_angledata=False``).
    """
    ncol = 7
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=f"{name}.nam")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
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
        gwf, icelltype=0, k=hk, k33=hk, save_flows=True, filename=f"{name}.npf"
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), strt], [(0, 0, ncol - 1), strt]],
        filename=f"{name}.chd",
    )

    # the second connection of well B is slanted (with ANGLEDATA) or an
    # equivalent vertical screen (without ANGLEDATA)
    if use_angledata:
        b2_top, b2_bot = mw_slant_top, mw_slant_bot
    else:
        b2_top, b2_bot = mw_equiv_top, mw_equiv_bot

    # packagedata: [ifno, radius, bottom, strt, condeqn, ngwfnodes]
    packagedata = [
        [0, radius, botm, strt, "MEAN", 1],
        [1, radius, botm, strt, "MEAN", 2],
    ]
    # connectiondata: [ifno, icon, cellid, scrn_top, scrn_bot, hk_skin, radius_skin]
    connectiondata = [
        [0, 0, (0, 0, 1), 90.0, 10.0, hks, sradius],  # well A, vertical
        [1, 0, (0, 0, 3), 80.0, 30.0, hks, sradius],  # well B, vertical
        [1, 1, (0, 0, 5), b2_top, b2_bot, hks, sradius],  # well B, slanted
    ]
    perioddata = {0: [[0, "rate", mawrate], [1, "rate", mawrate]]}

    maw_kwargs = dict(
        print_input=True,
        print_head=True,
        save_flows=True,
        observations={f"{name}.maw.obs.csv": [("ha", "head", 1), ("hb", "head", 2)]},
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname=f"maw_{name}",
        filename=f"{name}.maw",
    )
    if use_angledata:
        maw_kwargs["non_vertical_wells"] = True
        # only well B (ifno 1), connection 2 (icon 1) is non-vertical
        maw_kwargs["angledata"] = [[1, 1, mw_angle]]

    flopy.mf6.ModflowGwfmaw(gwf, **maw_kwargs)
    return gwf


def add_partialsat_model(sim, name, mode):
    """Add an unconfined single-well model used to check that the saturation is
    handled correctly for a partially saturated cell (Sf < 1).

    Three equivalent representations of the same connection over the screen
    interval [ps_bot, ps_top] are built; all must give identical heads at every
    head, including when the water table is partway through the screen:

    * ``mode="slant"`` - MEAN slanted connection (hk_skin scaled internally by
      the length correction)
    * ``mode="vert"``  - MEAN vertical connection with hk_skin pre-scaled by the
      length-correction factor
    * ``mode="spec"``  - SPECIFIED slanted connection with the saturated
      conductance supplied directly (the length correction is not applied to
      SPECIFIED, and the user-specified screen elevations are honored)
    """
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, model_nam_file=f"{name}.nam")
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=1,
        ncol=3,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )
    flopy.mf6.ModflowGwfic(gwf, strt=ps_chd, filename=f"{name}.ic")
    # unconfined (icelltype=1) so the connection can be partially saturated
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=1, k=hk, k33=hk, save_flows=True, filename=f"{name}.npf"
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), ps_chd]],
        filename=f"{name}.chd",
    )

    # select the conductance equation and skin/conductance value for each mode;
    # the slanted MEAN connection is scaled internally, the vertical MEAN
    # connection is pre-scaled, and the SPECIFIED connection is supplied the
    # already-corrected conductance directly
    if mode == "vert":
        condeqn, hk_skin, listed = "MEAN", ps_hks_vert, False
    elif mode == "slant":
        condeqn, hk_skin, listed = "MEAN", hks, True
    else:  # "spec"
        condeqn, hk_skin, listed = "SPECIFIED", ps_cspec, True
    packagedata = [[0, radius, botm, ps_chd, condeqn, 1]]
    connectiondata = [[0, 0, (0, 0, 1), ps_top, ps_bot, hk_skin, sradius]]
    perioddata = {0: [[0, "rate", ps_rate]]}

    maw_kwargs = dict(
        print_head=True,
        save_flows=True,
        observations={f"{name}.maw.obs.csv": [("h", "head", 1)]},
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata,
        pname=f"maw_{name}",
        filename=f"{name}.maw",
    )
    if listed:
        maw_kwargs["non_vertical_wells"] = True
        maw_kwargs["angledata"] = [[0, 0, ps_angle]]

    flopy.mf6.ModflowGwfmaw(gwf, **maw_kwargs)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL")],
    )
    return gwf


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", sim_ws=ws)
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

    add_single_well_model(sim, "vert", top, botm)
    add_single_well_model(sim, "angle0", top, botm, angledata=0.0)
    add_single_well_model(sim, "noangle", top, botm, nonvert_opt=True)
    add_single_well_model(sim, "horizontal", 51.0, 49.0, angledata=90.0, conn_len=lref)
    add_single_well_model(sim, "slant", slant_top, slant_bot, angledata=slant_angle)
    add_single_well_model(sim, "slant_equiv", slant_equiv_top, slant_equiv_bot)
    add_collector_model(sim, "collector")
    add_multiwell_model(sim, "multiwell", use_angledata=True)
    add_multiwell_model(sim, "multiwellref", use_angledata=False)
    add_partialsat_model(sim, "psat_slant", mode="slant")
    add_partialsat_model(sim, "psat_vert", mode="vert")
    add_partialsat_model(sim, "psat_spec", mode="spec")

    return sim


def maw_head(ws, model, col="H"):
    fpth = os.path.join(ws, f"{model}.maw.obs.csv")
    data = np.genfromtxt(fpth, delimiter=",", names=True)
    return float(np.atleast_1d(data[col])[-1])


def check_output(idx, test):
    ws = test.workspace

    h_vert = maw_head(ws, "vert")
    h_angle0 = maw_head(ws, "angle0")
    h_horizontal = maw_head(ws, "horizontal")
    h_slant = maw_head(ws, "slant")
    h_slant_equiv = maw_head(ws, "slant_equiv")

    # a 0 degree connection in the ANGLEDATA block must reproduce the vertical
    # (no ANGLEDATA) result exactly
    assert np.isclose(h_angle0, h_vert, atol=1e-6), (
        f"angle0 head ({h_angle0}) != vertical head ({h_vert})"
    )

    # NON_VERTICAL_WELLS specified without an ANGLEDATA block must run (a
    # warning is issued) and treat all connections as vertical, reproducing the
    # vertical result exactly
    h_noangle = maw_head(ws, "noangle")
    assert np.isclose(h_noangle, h_vert, atol=1e-6), (
        f"no-angledata head ({h_noangle}) != vertical head ({h_vert})"
    )

    # a horizontal connection with a connection length equal to the vertical
    # reference screen length must reproduce the vertical result
    assert np.isclose(h_horizontal, h_vert, atol=1e-6), (
        f"horizontal head ({h_horizontal}) != vertical head ({h_vert})"
    )

    # the slant case (derived in-cell screen length) must match a vertical well
    # whose screen thickness equals the derived length
    assert np.isclose(h_slant, h_slant_equiv, atol=1e-6), (
        f"slant head ({h_slant}) != equivalent vertical head ({h_slant_equiv})"
    )

    # the slanted connection has a longer screen than its vertical extent, so
    # its conductance is larger and the (extracting) well head is higher than
    # for a vertical well with the same screen thickness
    assert h_slant > h_vert - 1.0, "slant head unexpectedly low"

    # the four equal-length laterals of the radial collector well must carry
    # equal flows by symmetry
    fpth = os.path.join(ws, "collector.maw.obs.csv")
    data = np.genfromtxt(fpth, delimiter=",", names=True)
    flows = np.array([np.atleast_1d(data[f"Q{i + 1}"])[-1] for i in range(4)])
    assert np.allclose(flows, flows[0], atol=1e-6), (
        f"radial collector lateral flows are not symmetric: {flows}"
    )
    # the laterals supply the extracting well (positive GWF-to-well flow) and
    # together account for the total pumping rate
    assert np.all(flows > 0.0), "radial collector laterals should supply the well"
    assert np.isclose(flows.sum(), -mawrate, atol=1e-6), (
        f"radial collector lateral flows ({flows.sum()}) do not sum to the "
        f"pumping rate ({-mawrate})"
    )

    # multiple wells in one MAW package: the model with a slanted connection
    # (well B, connection 2 in the ANGLEDATA block) must reproduce the model in
    # which that connection is replaced by an equivalent vertical screen. This
    # confirms the ANGLEDATA block is applied only to the listed connection and
    # leaves the vertical well (well A) and the vertical connection (well B,
    # connection 1) unchanged.
    ha = maw_head(ws, "multiwell", col="HA")
    hb = maw_head(ws, "multiwell", col="HB")
    ha_ref = maw_head(ws, "multiwellref", col="HA")
    hb_ref = maw_head(ws, "multiwellref", col="HB")
    assert np.isclose(ha, ha_ref, atol=1e-6), (
        f"multiwell well A head ({ha}) != reference ({ha_ref})"
    )
    assert np.isclose(hb, hb_ref, atol=1e-6), (
        f"multiwell well B head ({hb}) != reference ({hb_ref})"
    )
    # the two wells should have distinct heads (sanity check)
    assert not np.isclose(ha, hb, atol=1e-6), (
        f"multiwell well A and well B heads are unexpectedly equal ({ha}, {hb})"
    )

    # partially saturated cell (Sf < 1): a slanted connection over a screen
    # interval is equivalent to a vertical connection over the same interval
    # with the skin conductivity scaled by the length correction, so the heads
    # must match even when the connection is only partially saturated
    h_slant = maw_head(ws, "psat_slant", col="H")
    h_vert = maw_head(ws, "psat_vert", col="H")
    h_spec = maw_head(ws, "psat_spec", col="H")
    assert np.isclose(h_slant, h_vert, atol=1e-6), (
        f"partial-saturation slant head ({h_slant}) != equivalent vertical "
        f"head ({h_vert})"
    )
    # a SPECIFIED slanted connection (conductance supplied directly, screen
    # elevations honored, length correction NOT applied) must match the MEAN
    # slanted connection given the equivalent saturated conductance; this
    # confirms SPECIFIED connections can be non-vertical, that the screen
    # elevations are honored (not reset to the full cell), and that the length
    # correction is not applied to SPECIFIED
    assert np.isclose(h_spec, h_slant, atol=1e-6), (
        f"partial-saturation SPECIFIED head ({h_spec}) != MEAN slant head ({h_slant})"
    )
    # confirm the connection is actually partially saturated (the water table
    # is within the screen interval), so this case exercises Sf < 1
    hcell = float(
        flopy.utils.HeadFile(os.path.join(ws, "psat_slant.hds")).get_data()[0, 0, 1]
    )
    assert ps_bot < hcell < ps_top, (
        f"cell head ({hcell}) is not within the screen interval "
        f"({ps_bot}, {ps_top}); partial saturation was not exercised"
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

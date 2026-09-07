"""
Test the ULTIMATE advection scheme where flow runs oblique to the grid.

A step front is advected across a square grid by a uniform flow field, once
along the grid and once at 45 degrees, at the same cell size, the same velocity
component along each axis and the same Courant number per direction.  The
exact answer is a step that has travelled a known distance, so both cases can be
measured directly.

The oblique case is the one that matters.  Along the grid every transverse and
twist term of the reconstruction vanishes identically and the scheme reduces to
one-dimensional QUICKEST, so the aligned case only confirms that the
multidimensional terms stay out of the way when they should.  At 45 degrees they
are as large as they get, and they are what holds the solution together: with
the face-normal terms alone this problem reaches a concentration of 2.8 and
minus 0.97 while still satisfying the Courant constraint, because the
one-dimensional limiter bounds each face on the assumption that the face carries
the only flux leaving the cell, and an obliquely flowing cell discharges through
two faces at once.  With the full reconstruction the same run stays inside
[0, 1.06] and its error halves.

Some overshoot survives, as it does in MT3DMS, which overshoots slightly more
than this on the same problem.  Multidimensional monotonicity is not guaranteed
by a limiter derived in one dimension, so the test bounds the overshoot rather
than forbidding it.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ult2d_align", "ult2d_obliq"]

ncell = 40  # cells per side
delta = 1.0
porosity = 0.25
hydraulic_conductivity = 1.0
gradient = 0.01  # head gradient along each grid axis
# seepage velocity component along each axis
velocity = hydraulic_conductivity * gradient / porosity
total_time = 375.0  # front advances velocity*total_time = 15 cells
courant = 0.45  # per direction, so the oblique case sums to 0.90


def oblique(name):
    return name.endswith("obliq")


def build_models(idx, name, test):
    dt = courant * delta / velocity
    # ATS_PERCEL bounds the Courant number summed over a cell's outgoing faces,
    # which is twice the per-direction value when flow is oblique
    percel = min(0.98, (2.0 if oblique(name) else 1.0) * courant + 0.03)

    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, nper=1, perioddata=[(total_time, round(total_time / dt), 1.0)]
    )
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, dt, dt * 1.0e-9, dt, 2.0, 2.0)],
        filename=f"{name}.ats",
    )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-8,
        linear_acceleration="CG",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=ncell,
        ncol=ncell,
        delr=delta,
        delc=delta,
        top=1.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=hydraulic_conductivity, save_specific_discharge=True
    )

    # prescribed heads on the whole perimeter reproduce a uniform gradient; the
    # inflow edges carry the inlet concentration as an auxiliary variable
    chd = []
    for i in range(ncell):
        for j in range(ncell):
            if not (i in (0, ncell - 1) or j in (0, ncell - 1)):
                continue
            x, y = (j + 0.5) * delta, (i + 0.5) * delta
            head = 10.0 - gradient * (x + y if oblique(name) else x)
            inflow = (j == 0 or i == 0) if oblique(name) else (j == 0)
            chd.append([(0, i, j), head, 1.0 if inflow else 0.0])
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: chd},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )

    gwtname = "gwt_" + name
    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-8,
        inner_dvclose=1.0e-8,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=1,
        nrow=ncell,
        ncol=ncell,
        delr=delta,
        delc=delta,
        top=1.0,
        botm=0.0,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="ultimate", ats_percel=percel)
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )
    return sim, None


def check_output(idx, name, test):
    gwtname = "gwt_" + name
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()[0]
    assert conc.shape == (ncell, ncell)

    xc = (np.arange(ncell) + 0.5) * delta
    x, y = np.meshgrid(xc, xc)
    front = velocity * total_time
    exact = np.where((np.minimum(x, y) if oblique(name) else x) < front, 1.0, 0.0)

    # the interior, so the prescribed boundary does not dominate the measure
    sel = (slice(2, ncell - 2), slice(2, ncell - 2))
    interior = conc[sel]

    # the solution must stay within the range of the data it was built from.
    # this is the assertion that fails without the multidimensional limiter:
    # with the face-normal terms alone the oblique case reaches 2.8 and -0.97
    cmin, cmax = interior.min(), interior.max()
    if oblique(name):
        # the one-dimensional limiter does not guarantee multidimensional
        # monotonicity, so the bound documents the overshoot rather than
        # forbidding it; measured 1.052 here, and MT3DMS overshoots more
        assert cmax < 1.06, f"oblique overshoot, max = {cmax}"
        assert cmin > -1.0e-5, f"oblique undershoot, min = {cmin}"
    else:
        # along the grid the transverse terms vanish and the limiter is exact
        assert cmax < 1.0 + 1.0e-5, f"aligned overshoot, max = {cmax}"
        assert cmin > -1.0e-5, f"aligned undershoot, min = {cmin}"
    print(f"{name}: min = {cmin:.6f}, max = {cmax:.6f}")

    # advection is explicit, so mass balance must be exact
    fpth = os.path.join(test.workspace, f"{gwtname}.lst")
    balance_error = None
    with open(fpth) as lst_file:
        for line in lst_file:
            if line.lstrip().startswith("PERCENT"):
                balance_error = float(line.split()[3])
    assert balance_error is not None, "balance error not found in listing file"
    assert abs(balance_error) < 1.0e-6, f"balance error = {balance_error}"

    error = np.abs(interior - exact[sel]).mean()
    if oblique(name):
        # measured 0.0325; upstream weighting gives 0.086 and the face-normal
        # terms on their own give 0.060
        assert error < 0.040, f"oblique front too diffuse, L1 per cell = {error}"
    else:
        # measured 0.0255 on this grid; upstream weighting gives 0.062
        assert error < 0.035, f"aligned front too diffuse, L1 per cell = {error}"


@pytest.mark.parametrize("name", cases)
def test_mf6model(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(0, name, t),
        check=lambda t: check_output(0, name, t),
    )
    test.run()

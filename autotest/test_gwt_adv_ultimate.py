"""
Test the explicit third-order ULTIMATE advection scheme on a one-dimensional
problem of pure advection through a row of square cells.  A step in
concentration is introduced at the inflow end and advected at a constant
velocity, so the exact answer is a step that has moved a known distance.

The same problem is set up three times, oriented along columns, along rows, and
along layers, to confirm that the scheme treats all three grid directions
identically.  Each orientation is run twice, at a Courant number of 1/4 and at a
Courant number of 1.  At a Courant number of 1 the QUICKEST reconstruction
underlying the scheme is exact, so the front must be resolved with no smearing
at all.

One case, "ult_dtadj0", sets DTADJ to zero, which is a legal input value
meaning "do not adjust the time step by the multiplier".  A package that submits
an absolute stability limit must still have that limit applied: ATS assigns the
submitted length to delt outright rather than scaling it, so DTADJ has no
bearing on it.  The case starts from a DT0 far below the stable step and asserts
that the step climbs to the stability limit anyway, which is what fails when the
submission is discarded.  The same path carries the SFR Package's kinematic wave
limit, so this guards more than the advection schemes.

One case, "ult_pt", puts a vertical pass-through cell (idomain = -1) in the
middle of the column.  Such a cell is removed from the solution but lets flow
cross it, so the cells either side are connected directly and are more than one
layer apart.  Finding the second upstream cell by arithmetic on layer indices
lands on the wrong cell there, and the scheme quietly falls back to first-order
upwinding around the gap.  MODFLOW 6 builds that connection from the two
neighbours' own thicknesses alone, so the column is geometrically identical to
one with no gap and must give an identical answer.

This case runs at a Courant number of a quarter rather than one.  At a Courant
number of one the update is a pure shift by one cell and first-order upwinding
is exact as well, so the fallback would be invisible.

Two further cases add linear sorption with a retardation factor of two.  The
scheme reconstructs the face value by tracing the front backward over the time
step, so it has to use the retarded velocity at which the front actually
advances rather than the water velocity.  If it does, these cases behave exactly
like the unretarded ones with the front at half the distance, and the Courant
number of one case is still exact.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "ult_col",
    "ult_row",
    "ult_lay",
    "ult_col_c1",
    "ult_row_c1",
    "ult_lay_c1",
    "ult_r2col",
    "ult_r2col_c1",
    "ult_dtadj0",
    "ult_pt",
]

# the layer index made a vertical pass-through cell in the "ult_pt" case
kpass = 40

# grid and flow are arranged so the seepage velocity is 10 length units per
# time unit and the front has advanced 50 cells at the end of the simulation
ncell = 100
delta = 1.0
porosity = 0.1
velocity = 1.0 / (porosity * delta * delta)
perlen = 5.0
front = velocity * perlen


def passthrough(name):
    """Case that puts a vertical pass-through cell in the middle of the column."""
    return "_pt" in name


def orientation(name):
    """Return grid shape and a cell index function for this case."""
    if "_col" in name:
        return (1, 1, ncell), lambda i: (0, 0, i)
    elif "_row" in name:
        return (1, ncell, 1), lambda i: (0, i, 0)
    elif passthrough(name):
        # one extra layer, which idomain removes from the solution, so the
        # active cells are still ncell.  Active cell i is layer i until the
        # pass-through and layer i+1 after it.
        return (ncell + 1, 1, 1), lambda i: (i if i < kpass else i + 1, 0, 0)
    else:
        return (ncell, 1, 1), lambda i: (i, 0, 0)


def idomain_kwargs(name):
    """DIS input marking the pass-through cell, if this case has one."""
    if not passthrough(name):
        return {}
    idomain = np.ones((ncell + 1, 1, 1), dtype=int)
    idomain[kpass, 0, 0] = -1
    return {"idomain": idomain}


def courant(name):
    return 1.0 if name.endswith("_c1") else 0.25


def dtadj0(name):
    """Case that checks a submitted stability limit survives DTADJ = 0."""
    return name.endswith("_dtadj0")


def retardation(name):
    """Retardation factor for this case, achieved with linear sorption."""
    return 2.0 if "_r2" in name else 1.0


def sorption_kwargs(name):
    """MST sorption input giving the requested retardation factor."""
    r = retardation(name)
    if r == 1.0:
        return {}
    bulk_density = 1.0
    return dict(
        sorption="linear",
        bulk_density=bulk_density,
        distcoef=(r - 1.0) * porosity / bulk_density,
    )


def build_models(idx, name, test):
    (nlay, nrow, ncol), cellid = orientation(name)
    # the stability limit is set by the retarded front velocity, so the time
    # step that gives the target Courant number grows with the retardation
    dt = courant(name) * retardation(name) * delta / velocity

    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=1,
        perioddata=[(perlen, round(perlen / dt), 1.0)],
    )
    # the scheme is explicit, so the time step must be under ATS control.
    # dtmax is pinned to dt so that every step is at the target Courant number.
    # dtadj must exceed one: ATS discards every submitted time step constraint
    # when dtadj is zero, which would leave the scheme with no stability limit.
    if dtadj0(name):
        # start far below the stable step, with no multiplier to grow it; only
        # the limit submitted by the ADV Package can raise it
        tdis.ats.initialize(
            maxats=1,
            perioddata=[(0, dt * 1.0e-3, dt * 1.0e-9, dt, 0.0, 0.0)],
            filename=f"{name}.ats",
        )
    else:
        tdis.ats.initialize(
            maxats=1,
            perioddata=[(0, dt, 1.0e-8, dt, 2.0, 2.0)],
            filename=f"{name}.ats",
        )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-6,
        inner_dvclose=1.0e-6,
        rcloserecord=1.0e-6,
        linear_acceleration="CG",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delta,
        delc=delta,
        top=nlay * delta,
        botm=[(nlay - k - 1) * delta for k in range(nlay)],
        **idomain_kwargs(name),
    )
    flopy.mf6.ModflowGwfic(gwf, strt=nlay * delta)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=1.0, k33=1.0, save_specific_discharge=True
    )
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[cellid(ncell - 1), 0.0]]})
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data={0: [[cellid(0), 1.0, 1.0]]},
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )

    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(sim, model_type="gwt6", modelname=gwtname)
    gwt.name_file.save_flows = True
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-6,
        inner_dvclose=1.0e-6,
        rcloserecord=1.0e-6,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delta,
        delc=delta,
        top=nlay * delta,
        botm=[(nlay - k - 1) * delta for k in range(nlay)],
        **idomain_kwargs(name),
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="ultimate", ats_percel=courant(name))
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, **sorption_kwargs(name))
    flopy.mf6.ModflowGwtssm(gwt, sources=[("WEL-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "LAST")],
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
    conc = cobj.get_data().flatten()
    if passthrough(name):
        # the pass-through cell is not in the solution and carries no value
        conc = conc[np.abs(conc) < 1.0e29]
    assert conc.shape == (ncell,)

    distance = (np.arange(ncell) + 0.5) * delta
    exact = np.where(distance < front / retardation(name), 1.0, 0.0)

    # the limiter must not allow the solution outside the range of the data
    assert conc.min() > -1.0e-8, f"undershoot: min concentration {conc.min()}"
    assert conc.max() < 1.0 + 1.0e-8, f"overshoot: max concentration {conc.max()}"

    # the profile must be monotonic, as there is no source of oscillation
    assert np.all(np.diff(conc) < 1.0e-8), "concentration profile is not monotonic"

    # advection is explicit, so mass balance must be exact
    fpth = os.path.join(test.workspace, f"{gwtname}.lst")
    balance_error = None
    with open(fpth) as lst_file:
        for line in lst_file:
            if line.lstrip().startswith("PERCENT"):
                balance_error = float(line.split()[3])
    assert balance_error is not None, "balance error not found in listing file"
    assert abs(balance_error) < 1.0e-6, f"balance error = {balance_error}"

    if dtadj0(name):
        # the run would take a thousand times more steps if the limit submitted
        # by the ADV Package were discarded because DTADJ is zero
        nstep = 0
        with open(os.path.join(test.workspace, "mfsim.lst")) as f:
            for line in f:
                if "ATS: time step set to" in line:
                    nstep += 1
        assert nstep < 3 * perlen / (courant(name) * delta / velocity), (
            f"time step never rose above DT0: {nstep} steps taken"
        )

    error = np.abs(conc - exact).sum() * delta
    if courant(name) == 1.0:
        # at a Courant number of one the scheme is exact
        assert error < 1.0e-8, f"front is not exact, L1 error = {error}"
    else:
        # otherwise the front is spread over a few cells; the upstream scheme
        # gives 6.30 and the TVD schemes give 3.14 on this problem
        assert error < 1.5, f"front is too diffuse, L1 error = {error}"
        assert error > 0.5, f"front is suspiciously sharp, L1 error = {error}"

    if passthrough(name):
        # the column is geometrically identical to one with no gap, so the
        # answer has to be identical too.  The equivalent ult_lay case gives
        # 1.3099; finding the second upstream cell by arithmetic on layer
        # indices misses across the gap, drops to first-order upwinding there
        # and gives 1.5964.
        assert error < 1.35, (
            f"the pass-through cell degraded the front, L1 error = {error}"
        )


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

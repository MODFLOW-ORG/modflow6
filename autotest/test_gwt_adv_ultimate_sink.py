"""
Test the ULTIMATE time step limit where flow converges on an extraction well.

A row of cells with equal constant heads at both ends and an extraction well in
the middle, so water flows toward the well from both sides. The well cell
receives water through two faces and loses none through its faces, so the
stable time step is set by the outflow of the other cells. Water entering at
the left end carries a concentration of one and at the right end zero, so the
well cell settles at one half.

Cases:
  - ult_sink     : no sorption; the stable time step is one day
  - ult_sink_ret : a retardation factor of two; the stable time step is two days
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ult_sink", "ult_sink_ret"]

# each half of the row carries half the well rate toward the well, giving a
# seepage velocity of one cell per day and an unretarded Courant number of one
# at a one-day time step
ncol = 41
iwell = 20
delta = 1.0
porosity = 0.25
rate = 0.5
perlen = 80.0
dt0 = 0.01


def retardation(name):
    """Return the retardation factor of a case."""
    return 2.0 if name.endswith("_ret") else 1.0


def stable_dt(name):
    """Return the time step at a Courant number of one."""
    pore_volume = delta * delta * delta * porosity * retardation(name)
    return pore_volume / (0.5 * rate)


def build_models(idx, name, test):
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    tdis = flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, 1, 1.0)]
    )
    # a short first step with no multiplier, so the limit submitted by the ADV
    # Package sets every later step
    tdis.ats.initialize(
        maxats=1,
        perioddata=[(0, dt0, 1.0e-9, perlen, 0.0, 0.0)],
        filename=f"{name}.ats",
    )
    dis = dict(nlay=1, nrow=1, ncol=ncol, delr=delta, delc=delta, top=delta, botm=0.0)

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-10,
        inner_dvclose=1.0e-10,
        rcloserecord=1.0e-10,
        linear_acceleration="CG",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])
    flopy.mf6.ModflowGwfdis(gwf, **dis)
    flopy.mf6.ModflowGwfic(gwf, strt=10.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), 10.0, 1.0], [(0, 0, ncol - 1), 10.0, 0.0]]},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfwel(gwf, stress_period_data={0: [[(0, 0, iwell), -rate]]})

    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(sim, model_type="gwt6", modelname=gwtname)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1.0e-10,
        inner_dvclose=1.0e-10,
        rcloserecord=1.0e-10,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(gwt, **dis)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="ultimate", ats_percel=1.0)
    sorption = {}
    if retardation(name) != 1.0:
        sorption = dict(
            sorption="linear",
            bulk_density=1.0,
            distcoef=(retardation(name) - 1.0) * porosity,
        )
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, **sorption)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "LAST")],
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
    assert conc.shape == (ncol,)

    # every step after the first is set by the ADV Package; counting the well
    # cell inflow, which leaves through the well rather than a cell face, would
    # halve the step
    submitted = []
    nstep = 0
    with open(os.path.join(test.workspace, "mfsim.lst")) as f:
        for line in f:
            if "ATS: time step set to" in line:
                nstep += 1
            if "submitted a preferred time step" in line and "ADV" in line:
                submitted.append(float(line.split()[-1]))
    assert submitted, "ADV Package submitted no time step"
    dt = stable_dt(name)
    assert np.allclose(submitted, dt, rtol=1.0e-10), (
        f"ADV Package submitted {min(submitted)} to {max(submitted)}, expected {dt}"
    )
    expected = 1 + int(np.ceil((perlen - dt0) / dt - 1.0e-9))
    assert nstep == expected, f"took {nstep} time steps, expected {expected}"

    # the front has reached the well, so the inflow sides are uniform and the
    # well cell is the flow-weighted mean of the two
    assert conc.min() > -1.0e-10, f"undershoot: min concentration {conc.min()}"
    assert conc.max() < 1.0 + 1.0e-10, f"overshoot: max concentration {conc.max()}"
    assert np.allclose(conc[:iwell], 1.0, atol=1.0e-8), "left side is not one"
    assert np.allclose(conc[iwell + 1 :], 0.0, atol=1.0e-8), "right side is not zero"
    assert abs(conc[iwell] - 0.5) < 1.0e-8, f"well cell is {conc[iwell]}"

    # advection is explicit, so mass balance must be exact
    balance_error = None
    with open(os.path.join(test.workspace, f"{gwtname}.lst")) as f:
        for line in f:
            if line.lstrip().startswith("PERCENT"):
                balance_error = float(line.split()[3])
    assert balance_error is not None, "balance error not found in listing file"
    assert abs(balance_error) < 1.0e-6, f"balance error = {balance_error}"


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

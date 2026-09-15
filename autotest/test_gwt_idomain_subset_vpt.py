"""
Tests a GWT model that excludes a vertical pass-through cell of the GWF model.

A three-layer flow model has a vertical pass-through cell (IDOMAIN of -1) in the
middle layer, so the cells above and below it are connected.  The transport
model gives that cell an IDOMAIN of 0, so both models have the same number of
active cells but the transport model does not connect the cells above and below.

Cases:
  - idomvpt : flow through the pass-through connection is carried by the flow
              imbalance correction; the transport budget must close and
              concentrations must stay within the source range.
"""

import os
import re

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["idomvpt"]

nlay, nrow, ncol = 3, 1, 10
delr, delc = 1.0, 1.0
top, botm = 3.0, [2.0, 1.0, 0.0]
hk = 1.0
porosity = 0.1
perlen, nstp = 10.0, 50

# middle-layer column that is a pass-through cell in the flow model
icolvpt = 4


def build_models(idx, test):
    name = cases[idx]
    nouter, ninner = 100, 300
    hclose, rclose = 1e-8, 1e-8

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=1, perioddata=[(perlen, nstp, 1.0)]
    )

    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(ims, [gwf.name])
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[1, 0, icolvpt] = -1
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(
        gwf, icelltype=0, k=hk, k33=hk, save_specific_discharge=True
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), 1.0, 1.0], [(2, 0, ncol - 1), 0.0, 0.0]]},
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(sim, model_type="gwt6", modelname=gwtname)
    gwt.name_file.save_flows = True
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(ims, [gwt.name])
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[1, 0, icolvpt] = 0
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtadv(gwt, scheme="upstream")
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)
    flopy.mf6.ModflowGwtssm(gwt, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("BUDGET", "ALL")],
    )

    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    name = cases[idx]
    gwtname = "gwt_" + name

    # the transport budget must close
    with open(os.path.join(test.workspace, f"{gwtname}.lst")) as f:
        text = f.read()
    residual = np.array(re.findall(r"IN - OUT =\s*(\S+)", text), dtype=float)
    total = np.array(re.findall(r"TOTAL IN =\s*(\S+)", text), dtype=float)
    assert residual.size > 0, f"no budget written to {gwtname}.lst"
    pd = np.abs(residual) / np.where(total == 0.0, 1.0, total)
    assert pd.max() < 1e-6, (
        f"transport budget does not close, max discrepancy {pd.max()}"
    )

    # concentrations in active cells must stay within the source range
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()
    active = conc != 1e30
    assert np.all(conc[active] >= -1e-9), "negative concentration"
    assert np.all(conc[active] <= 1.0 + 1e-9), "concentration above the source"


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

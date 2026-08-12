"""
Tests MWT budget output when every well is inactive in the first stress period.

The GWF node number written to the advanced package budget was only assigned for
active features, so an inactive feature wrote whatever the uninitialized local
held. An idomain hole makes the node numbering reduced, so a stale node number is
converted through the user-node lookup on output.

Cases:
  - mwt_inactive : all wells INACTIVE in period 1 and ACTIVE afterwards; the GWF
                   budget entries must carry the connected cells in every period,
                   with zero flow while the wells are inactive.
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["mwt_inactive"]

nlay, nrow, ncol, nper = 3, 1, 5, 3
well_cols = (1, 3)

# user node numbers of the well connections, which are what the budget file holds
expected_nodes = sorted(k * nrow * ncol + j + 1 for j in well_cols for k in range(nlay))


def build_models(idx, test):
    name = cases[idx]
    gwfname, gwtname = "gwf_" + name, "gwt_" + name
    botm = [-1.0, -2.0, -3.0]

    # the idomain hole gives the model a reduced node numbering
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    idomain[nlay - 1, 0, 2] = 0

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(10.0, 5, 1.0)] * nper
    )

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-8,
        inner_dvclose=1e-9,
        linear_acceleration="bicgstab",
        filename=f"{gwfname}.ims",
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)
    sim.register_ims_package(imsgwf, [gwfname])
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=0.0,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=1e-5, sy=0.1, transient={0: True})
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[
            [(0, 0, 0), 0.0, 0.0],
            [(0, 0, ncol - 1), 0.0, 0.0],
        ],
        auxiliary="CONCENTRATION",
        pname="CHD-1",
    )

    packagedata, conndata, perdata = [], [], []
    for iw, jcol in enumerate(well_cols):
        packagedata.append([iw, 0.1, -3.0, 0.0, "THIEM", nlay])
        for k in range(nlay):
            conndata.append(
                [
                    iw,
                    k,
                    (k, 0, jcol),
                    0.0 if k == 0 else -float(k),
                    -float(k + 1),
                    1.0,
                    0.1,
                ]
            )
        perdata.append([iw, "rate", -0.5])
    flopy.mf6.ModflowGwfmaw(
        gwf,
        save_flows=True,
        packagedata=packagedata,
        connectiondata=conndata,
        perioddata=perdata,
        pname="MAW-1",
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname, save_flows=True)
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        complexity="moderate",
        outer_dvclose=1e-8,
        inner_dvclose=1e-9,
        linear_acceleration="bicgstab",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwtname])
    flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=0.0,
        botm=botm,
        idomain=idomain,
    )
    flopy.mf6.ModflowGwtic(gwt, strt=100.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtdsp(gwt, alh=0.1, ath1=0.01)
    flopy.mf6.ModflowGwtmwt(
        gwt,
        save_flows=True,
        budget_filerecord=f"{gwtname}.mwt.bud",
        concentration_filerecord=f"{gwtname}.mwt.bin",
        packagedata=[(iw, 0.0) for iw in range(len(well_cols))],
        mwtperioddata={
            0: [(iw, "STATUS", "INACTIVE") for iw in range(len(well_cols))],
            1: [(iw, "STATUS", "ACTIVE") for iw in range(len(well_cols))],
        },
        pname="MAW-1",
    )
    flopy.mf6.ModflowGwtssm(gwt, sources=[["CHD-1", "AUX", "CONCENTRATION"]])
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
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
    gwtname = "gwt_" + cases[idx]
    fpth = os.path.join(test.workspace, f"{gwtname}.mwt.bud")
    assert os.path.isfile(fpth)

    cbc = flopy.utils.CellBudgetFile(fpth, precision="double")
    for kstpkper in cbc.get_kstpkper():
        rec = cbc.get_data(kstpkper=kstpkper, text="GWF")[0]
        nodes = sorted(int(n) for n in rec["node2"])
        assert nodes == expected_nodes, (
            f"GWF budget nodes {nodes} at {kstpkper} do not match the well "
            f"connections {expected_nodes}"
        )
        # the wells are inactive for the whole of the first stress period
        if kstpkper[1] == 0:
            assert np.allclose(rec["q"], 0.0), (
                f"nonzero MWT-GWF flow {rec['q']} at {kstpkper} while inactive"
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

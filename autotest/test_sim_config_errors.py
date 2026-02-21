"""
Test to verify simulation configuration error checking.

Tests the following error conditions:
1. IMS6 solver paired with PRT model (should require EMS6)
2. EMS6 solver paired with GWF model (should require IMS6)
"""

import flopy
import pytest
from framework import TestFramework

simname = "cfgerr"
cases = [
    f"{simname}_ims_prt",
    f"{simname}_ems_gwf",
]


def build_ims_with_prt_sim(test):
    name = test.name
    ws = test.workspace
    mf6 = test.targets["mf6"]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    flopy.mf6.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=1,
        perioddata=[(1.0, 1, 1.0)],
    )

    prt_name = "prt1"
    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
        prt,
        pname="dis",
        nlay=1,
        nrow=10,
        ncol=10,
    )

    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=0.1)

    releasepts = [[0, 0, 0, 0, 0.5, 0.5, 0.5]]
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(releasepts),
        packagedata=releasepts,
        perioddata={0: [("FIRST",)]},
    )

    flopy.mf6.ModflowPrtoc(prt, pname="oc")
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFHEAD", "dummy.hds"),
            ("GWFBUDGET", "dummy.bud"),
        ],
    )

    # Add IMS solver (wrong - PRT needs EMS)
    ims = flopy.mf6.ModflowIms(sim, pname="ims")
    sim.register_solution_package(ims, [prt.name])

    return sim


def build_ems_with_gwf_sim(test):
    name = test.name
    ws = test.workspace
    mf6 = test.targets["mf6"]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    flopy.mf6.ModflowTdis(
        sim,
        pname="tdis",
        time_units="DAYS",
        nper=1,
        perioddata=[(1.0, 1, 1.0)],
    )

    gwf_name = "gwf1"
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name)

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=1,
        nrow=10,
        ncol=10,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True)

    chd_data = {0: [[(0, 0, 0), 1.0], [(0, 9, 9), 0.0]]}
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_data)

    # Add EMS solver (wrong - GWF needs IMS)
    ems = flopy.mf6.ModflowEms(sim, pname="ems", filename=f"{gwf_name}.ems")
    sim.register_solution_package(ems, [gwf.name])

    return sim


def build_models(test):
    if "ims_prt" in test.name:
        return build_ims_with_prt_sim(test)
    elif "ems_gwf" in test.name:
        return build_ems_with_gwf_sim(test)


def check_output(test):
    name = test.name
    buff = test.buffs[0] if len(test.buffs) > 0 else []

    if "ims_prt" in name:
        assert any("explicit model" in line.lower() and "ims6" in line.lower()
                   for line in buff)
    elif "ems_gwf" in name:
        assert any("numerical model" in line.lower() and "ems6" in line.lower()
                   for line in buff)


@pytest.mark.parametrize("name", cases)
def test_sln_config_errors(name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(t),
        check=lambda t: check_output(t),
        targets=targets,
        compare=None,
        xfail=True,
    )
    test.run()

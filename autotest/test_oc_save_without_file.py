"""
Verify that MODFLOW 6 terminates with an error when a model's Output
Control PERIOD block requests ``SAVE <rtype>`` but the corresponding
binary output file is not specified in the OC OPTIONS block.

Every output type managed by output control is exercised:

* GWF -- HEAD, BUDGET
* GWT -- CONCENTRATION, BUDGET
* GWE -- TEMPERATURE, BUDGET
* CHF -- STAGE, BUDGET

OLF shares the surface-water OC code path with CHF, so it is not tested
separately.
"""

import flopy
import pytest
from framework import TestFramework

cases = [
    ("gwf", "BUDGET"),
    ("gwf", "HEAD"),
    ("gwt", "BUDGET"),
    ("gwt", "CONCENTRATION"),
    ("gwe", "BUDGET"),
    ("gwe", "TEMPERATURE"),
    ("chf", "BUDGET"),
    ("chf", "STAGE"),
]


def _params():
    params = []
    for family, rtype in cases:
        # CHF/OLF models only build with a develop-mode binary
        marks = [pytest.mark.developmode] if family in ("chf", "olf") else []
        params.append(
            pytest.param(family, rtype, id=f"{family}-{rtype.lower()}", marks=marks)
        )
    return params


def build_gwf(ws, name, rtype):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple")
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(gwf, nlay=1, nrow=1, ncol=2)
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 1.0]])
    # no budget_filerecord / head_filerecord on purpose
    flopy.mf6.ModflowGwfoc(gwf, saverecord=[(rtype, "ALL")])
    return sim


def build_gwt(ws, name, rtype):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple")
    gwt = flopy.mf6.ModflowGwt(sim, modelname=name)
    flopy.mf6.ModflowGwtdis(gwt, nlay=1, nrow=1, ncol=2)
    flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)
    # no budget_filerecord / concentration_filerecord on purpose
    flopy.mf6.ModflowGwtoc(gwt, saverecord=[(rtype, "ALL")])
    return sim


def build_gwe(ws, name, rtype):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, complexity="simple")
    gwe = flopy.mf6.ModflowGwe(sim, modelname=name)
    flopy.mf6.ModflowGwedis(gwe, nlay=1, nrow=1, ncol=2)
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.2,
        density_water=1000.0,
        density_solid=2500.0,
        heat_capacity_water=4000.0,
        heat_capacity_solid=1000.0,
    )
    # no budget_filerecord / temperature_filerecord on purpose
    flopy.mf6.ModflowGweoc(gwe, saverecord=[(rtype, "ALL")])
    return sim


def build_chf(ws, name, rtype):
    sim = flopy.mf6.MFSimulation(sim_name=name, exe_name="mf6", sim_ws=ws)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    flopy.mf6.ModflowIms(sim, linear_acceleration="BICGSTAB")
    chf = flopy.mf6.ModflowChf(sim, modelname=name, save_flows=True)
    nreach, dx = 3, 1000.0
    vertices = [[j, j * dx, 0.0] for j in range(nreach + 1)]
    cell1d = [[j, 0.5, 2, j, j + 1] for j in range(nreach)]
    flopy.mf6.ModflowChfdisv1D(
        chf,
        nodes=nreach,
        nvert=len(vertices),
        width=50.0,
        bottom=0.0,
        idomain=1,
        vertices=vertices,
        cell1d=cell1d,
    )
    flopy.mf6.ModflowChfdfw(chf, manningsn=0.035, idcxs=0)
    flopy.mf6.ModflowChfsto(chf)
    flopy.mf6.ModflowChfic(chf, strt=1.0)
    xfraction = [0.0, 0.0, 1.0, 1.0]
    height = [100.0, 0.0, 0.0, 100.0]
    mannfraction = [1.0, 1.0, 1.0, 1.0]
    cxsdata = list(zip([0] * 4, xfraction, height, mannfraction))
    flopy.mf6.ModflowChfcxs(
        chf,
        nsections=1,
        npoints=4,
        packagedata=[(0, 4)],
        crosssectiondata=cxsdata,
    )
    # no budget_filerecord / stage_filerecord on purpose
    flopy.mf6.ModflowChfoc(chf, saverecord=[(rtype, "ALL")])
    return sim


BUILDERS = {
    "gwf": build_gwf,
    "gwt": build_gwt,
    "gwe": build_gwe,
    "chf": build_chf,
}


def build_models(test, family, rtype):
    return BUILDERS[family](test.workspace, f"oc{family}", rtype)


def check_output(test, rtype):
    with open(test.workspace / "mfsim.lst") as f:
        content = " ".join(f.read().split())

    expected = (
        f"REQUESTING TO SAVE {rtype} BUT {rtype} SAVE FILE NOT SPECIFIED. "
        f"{rtype} SAVE FILE MUST BE SPECIFIED IN OUTPUT CONTROL OPTIONS."
    )
    assert expected in content, f"missing expected OC error for {rtype}"


@pytest.mark.parametrize("family, rtype", _params())
def test_oc_save_without_file(family, rtype, function_tmpdir, targets):
    test = TestFramework(
        name=f"oc-savenofile-{family}-{rtype.lower()}",
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(t, family, rtype),
        check=lambda t: check_output(t, rtype),
        xfail=True,
    )
    test.run()

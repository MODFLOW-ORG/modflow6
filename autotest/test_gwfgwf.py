"""
Test a constant head boundary assigned to one model as
also functioning as a constant head boundary assigned to
a connected model.  The constant head term for model 1
should be calculated to include the flow from model 2.
Also, the flowja budget term for the constant head cell
should have a correct residual of zero in the diagonal
position.

   1  1 -1  gwf1
   -  -  -
   1  1  1  gwf2

We assert equality on qresidual being less than tolerance and
also that total constant head flow is correct.
"""

import math
import pathlib as pl

import flopy
import pytest
from compare import Comparison
from framework import TestFramework

cases = ["gwfgwf01", "gwfgwf01ifmod"]
ifmod = [False, True]


def build_models(idx, test):
    sim = get_sim(idx, test.workspace)
    return sim, None


def get_sim(idx, dir):
    name = cases[idx]

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1.0e-8, 1e-8, 0.97

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=dir
    )

    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=1,
        perioddata=[(1.0, 1, 1)],
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
    )

    gwf0 = add_model(sim, "gwf1", top=0.0, add_chd=True, add_rch=True)
    gwf1 = add_model(sim, "gwf2", top=-1.0, add_chd=False, add_rch=False)
    gwfgwf = add_gwfexchange(sim, idx)

    return sim


def add_model(sim, modelname, top, add_chd, add_rch):
    nlay, nrow, ncol = 1, 1, 3
    botm = top - 1.0
    gwf = flopy.mf6.ModflowGwf(sim, modelname=modelname, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=1.0,
        delc=1.0,
        top=top,
        botm=botm,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        save_flows=True,
        icelltype=0,
        k=1.0,
    )

    if add_chd:
        chdlist = [(0, 0, ncol - 1, 0.0)]
        chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdlist)

    if add_rch:
        rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.001)

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{modelname}.hds",
        budget_filerecord=f"{modelname}.cbc",
        budgetcsv_filerecord=f"{modelname}.cbc.csv",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    return gwf


def add_gwfexchange(sim, idx):
    ncol = 3
    delr = 1.0
    delc = 1.0
    dz = 1.0
    angldegx = 0.0
    cdist = 1.0
    gwfgwf_data = [
        [
            (0, 0, icol),
            (0, 0, icol),
            0,
            dz / 2.0,
            dz / 2.0,
            delr * delc,
            angldegx,
            cdist,
        ]
        for icol in range(ncol)
    ]
    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        save_flows=True,
        print_flows=True,
        nexg=len(gwfgwf_data),
        exgmnamea="gwf1",
        exgmnameb="gwf2",
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
        dev_interfacemodel_on=ifmod[idx],
    )
    return gwfgwf


def check_output(idx, test):
    sim = test.sims[0]
    check_model(sim, 0)
    check_model(sim, 1)


def check_model(sim, model_number):
    print(f"Checking model output (gwf{model_number})")
    gwf = sim.gwf[model_number]
    sim_ws = pl.Path(sim.sim_path)
    fpth = sim_ws / f"gwf{model_number + 1}.dis.grb"
    grb = flopy.mf6.utils.MfGrdFile(fpth)
    ia = grb.ia

    bobj = gwf.output.budget()
    # print(bobj.list_records())
    budget_records = bobj.get_data(kstpkper=(0, 0))
    nrecords = len(budget_records)
    print(f"detected {nrecords} items in budget file.")
    for idx in range(nrecords):
        budget_record = bobj.get_data(idx=idx)[0]
        print(budget_record)

    model_budget = gwf.output.budgetcsv().data
    pd = model_budget["PERCENT_DIFFERENCE"][0]
    print("percent difference: ", model_budget["PERCENT_DIFFERENCE"])
    errmsg = "Model percent difference is too large (pd)"
    assert pd < 1.0e-6, errmsg

    # check residual budget term in flowja diagonal position
    fja = bobj.get_data(text="FLOW-JA-FACE")[0].flatten()
    success = True
    atol = 1.0e-7
    for ipos in ia[:-1]:
        print(ipos, fja[ipos])
        if fja[ipos] > atol:
            success = False
    assert success, f"flowja residual larger than tolerance ({atol})"

    # ensure that the constant head outflow is equal to the
    # specified recharge
    if model_number == 0:
        chdflows = bobj.get_data(text="chd")[0]
        vchd = -chdflows["q"][0]
        vexgchd = model_budget["FLOW-JA-FACE-CHD(GWF-GWF_1)_OUT"][0]
        v = vchd + vexgchd
        errmsg = (
            f"Constant head outflow ({v}) is not equal to the "
            f"specified recharge inflow (0.002)."
        )
        assert math.isclose(v, 0.002), errmsg

    return


@pytest.mark.parametrize(
    "compare", ["auto", pytest.param(Comparison.FP4_ASCII, marks=pytest.mark.fp4)]
)
@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets, compare):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=compare,
    )
    test.run()

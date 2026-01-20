"""
Test for  GWFGWF exchange approach, with horizontal _and_ vertical
connections.  Ensure flowja(idiag) is set right for horizontal
and vertical connections.

layer 1:

        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 0 0 0 1 1
(H=1.0) 1 1 0 0 0 1 1 (H=0.0)   with the 9x9 refined grid inside
        1 1 0 0 0 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1

layer 2/3:

        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1
        1 1 1 1 1 1 1


There is no xt3d here or anything fancy. We just want to
ensure that vertical flows and horizontal flows are added
correctly to each model flowja diagonal terms.  This diagonal
term contains the flow residual for the cell.
"""

import os

import flopy
import numpy as np
import pytest
from conftest import project_root_path
from flopy.utils.lgrutil import Lgr
from framework import TestFramework

data_path = project_root_path / "autotest" / "data"
cases = [
    "gwfgwf_lgr_classic",
    "gwfgwf_lgr_ifmod",
    "gwfgwf_lgr_classic_binary",
    "gwfgwf_lgr_ifmod_binary",
]
ifmod = [False, True, False, True]
parent_name = "parent"
child_name = "child"
h_left = 1.0
h_right = 0.0
delr = 100.0
delc = 100.0
k11 = 1.0
k33 = 1.0


def get_model(idx, test):
    global child_domain
    global hclose

    name = cases[idx]

    # tdis period data
    nper = 1
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((1.0, 1, 1))

    # solver data
    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-12, 1e-3, 0.97
    h_start = 1.0

    # dis
    nlay, nrow, ncol = 3, 7, 7
    nlayc = 1

    row_s, row_e = 3, 5
    col_s, col_e = 3, 5

    ref_fct = 3
    nrowc = ref_fct * ((row_e - row_s) + 1)
    ncolc = ref_fct * ((col_e - col_s) + 1)

    idomain = np.ones((nlay, nrow, ncol))
    idomain[0, row_s - 1 : row_e, col_s - 1 : col_e] = 0

    delrc = delr / ref_fct
    delcc = delc / ref_fct
    tops = [0.0, -100.0, -200.0, -300.0]

    xoriginc = 2 * delr
    yoriginc = 2 * delc

    # boundary stress period data
    left_chd = [
        [(ilay, irow, 0), h_left]
        for ilay in range(1)  # apply chd only to top layer to drive vertical flow
        for irow in range(nrow)
    ]
    right_chd = [
        [(ilay, irow, ncol - 1), h_right]
        for ilay in range(1)  # apply chd only to top layer to drive vertical flow
        for irow in range(nrow)
    ]
    chd_data = left_chd + right_chd
    chd_spd = {0: chd_data}

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        memory_print_option="ALL",
    )

    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

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

    # The parent model:
    gwf = flopy.mf6.ModflowGwf(sim, modelname=parent_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=tops[0],
        botm=tops[1:],
        idomain=idomain,
    )
    ic = flopy.mf6.ModflowGwfic(gwf, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        xt3doptions=False,
        save_flows=True,
        icelltype=0,
        k=k11,
        k33=k33,
    )
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{parent_name}.hds",
        budget_filerecord=f"{parent_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # The child model:
    gwfc = flopy.mf6.ModflowGwf(sim, modelname=child_name, save_flows=True)
    dis = flopy.mf6.ModflowGwfdis(
        gwfc,
        nlay=nlayc,
        nrow=nrowc,
        ncol=ncolc,
        delr=delrc,
        delc=delcc,
        top=tops[0],
        botm=tops[1],
        xorigin=xoriginc,
        yorigin=yoriginc,
    )
    ic = flopy.mf6.ModflowGwfic(gwfc, strt=h_start)
    npf = flopy.mf6.ModflowGwfnpf(
        gwfc,
        save_specific_discharge=True,
        xt3doptions=False,
        save_flows=True,
        icelltype=0,
        k=k11,
        k33=k33,
    )
    oc = flopy.mf6.ModflowGwfoc(
        gwfc,
        head_filerecord=f"{child_name}.hds",
        budget_filerecord=f"{child_name}.cbc",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # LGR:
    nrowp = gwf.dis.nrow.get_data()
    ncolp = gwf.dis.ncol.get_data()
    delrp = gwf.dis.delr.array
    delcp = gwf.dis.delc.array
    topp = gwf.dis.top.array
    botmp = gwf.dis.botm.array
    idomainp = gwf.dis.idomain.array

    lgr = Lgr(
        3,
        nrowp,
        ncolp,
        delrp,
        delcp,
        topp,
        botmp,
        idomainp,
        ncpp=ref_fct,
        ncppl=[1, 0, 0],
    )

    exgdata = lgr.get_exchange_data(angldegx=True, cdist=True)

    if name.endswith("_binary"):
        exchangedata = {
            "factor": 1.0,
            "filename": "exg.bin",
            "data": exgdata,
            "binary": True,
        }
    else:
        exchangedata = exgdata

    gwfgwf = flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(exgdata),
        exgmnamea=parent_name,
        exgmnameb=child_name,
        exchangedata=exchangedata,
        print_flows=True,
        auxiliary=["ANGLDEGX", "CDIST"],
        dev_interfacemodel_on=ifmod[idx],
    )

    return sim


def build_models(idx, exdir):
    sim = get_model(idx, exdir)
    return sim, None


def check_output(idx, test):
    print("comparing heads  for child model to analytical result...")

    fpth = os.path.join(test.workspace, f"{child_name}.hds")
    hds_c = flopy.utils.HeadFile(fpth)
    heads_c = hds_c.get_data()

    fpth = os.path.join(test.workspace, f"{child_name}.dis.grb")
    grb_c = flopy.mf6.utils.MfGrdFile(fpth)

    # check flowja residual
    for mname in [parent_name, child_name]:
        print(f"Checking flowja residual for model {mname}")

        fpth = os.path.join(test.workspace, f"{mname}.dis.grb")
        grb = flopy.mf6.utils.MfGrdFile(fpth)
        ia = grb._datadict["IA"] - 1

        fpth = os.path.join(test.workspace, f"{mname}.cbc")
        assert os.path.isfile(fpth)
        cbb = flopy.utils.CellBudgetFile(fpth, precision="double")
        flow_ja_face = cbb.get_data(idx=0)
        assert len(flow_ja_face) > 0, (
            "Could not check residuals as flow-ja-face could not be found"
        )

        for fjf in flow_ja_face:
            fjf = fjf.flatten()
            res = fjf[ia[:-1]]
            errmsg = f"min or max residual too large {res.min()} {res.max()}"
            assert np.allclose(res, 0.0, atol=1.0e-6), errmsg


@pytest.mark.parametrize("idx, name", enumerate(cases))
@pytest.mark.developmode
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

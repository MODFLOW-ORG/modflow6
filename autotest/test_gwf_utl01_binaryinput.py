"""
Test reading of binary initial heads (float) and also binary icelltype (int).
1. Have binary data in a separate record for each layer
2. Have binary data in a single record for all layers
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["binary01", "binary02"]


def build_models(idx, test):
    nlay, nrow, ncol = 5, 6, 7
    nper = 1
    perlen = 1.0
    nstp = 1
    tsmult = 1.0
    steady = [True]
    lenx = 300.0
    delr = delc = lenx / float(nrow)
    botm = np.linspace(-1.0, -5.0, 5)
    hnoflo = 1e30
    hdry = -1e30
    hk = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-3, 1.0

    tdis_rc = []
    for _ in range(nper):
        tdis_rc.append((perlen, nstp, tsmult))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    # discretization
    # write top to a binary file
    text = "TOP"
    fname = "top.bin"
    pth = os.path.join(test.workspace, fname)
    f = open(pth, "wb")
    header = flopy.utils.BinaryHeader.create(
        bintype="HEAD",
        precision="double",
        text=text,
        nrow=nrow,
        ncol=ncol,
        ilay=1,
        pertim=1.0,
        totim=1.0,
        kstp=1,
        kper=1,
    )
    flopy.utils.Util2d.write_bin(
        (nrow, ncol),
        f,
        np.zeros((nrow, ncol), dtype=np.float64),
        header_data=header,
    )
    f.close()
    top = {
        "factor": 1.0,
        "filename": fname,
        "data": None,
        "binary": True,
        "iprn": 1,
    }

    # write botarr to binary file
    if idx == 0:
        botarr = []
        for k in range(nlay):
            text = f"BOTM_L{k + 1}"
            fname = f"botm.l{k + 1:02d}.bin"
            pth = os.path.join(test.workspace, fname)
            f = open(pth, "wb")
            header = flopy.utils.BinaryHeader.create(
                bintype="HEAD",
                precision="double",
                text=text,
                nrow=nrow,
                ncol=ncol,
                ilay=k + 1,
                pertim=1.0,
                totim=1.0,
                kstp=1,
                kper=1,
            )
            flopy.utils.Util2d.write_bin(
                (nrow, ncol),
                f,
                np.ones((nrow, ncol), dtype=np.float64) * botm[k],
                header_data=header,
            )
            f.close()
            botarr.append(
                {
                    "factor": 1.0,
                    "filename": fname,
                    "data": None,
                    "binary": True,
                    "iprn": 1,
                }
            )
    elif idx == 1:
        fname = "botm.bin"
        pth = os.path.join(test.workspace, fname)
        f = open(pth, "wb")
        tarr = np.ones((nlay, nrow, ncol), dtype=np.float64)
        for k in range(nlay):
            tarr[k, :, :] = botm[k]
        tarr = tarr.flatten()
        header = flopy.utils.BinaryHeader.create(
            bintype="HEAD",
            precision="double",
            text="BOTM",
            nrow=1,
            ncol=ncol * nrow * nlay,
            ilay=1,
            pertim=1.0,
            totim=1.0,
            kstp=1,
            kper=1,
        )
        flopy.utils.Util2d.write_bin((nrow, ncol), f, tarr, header_data=header)
        f.close()
        botarr = {
            "factor": 1.0,
            "filename": fname,
            "data": None,
            "binary": True,
            "iprn": 1,
        }

    # write idomain to binary file
    if idx == 0:
        idomain = []
        for k in range(nlay):
            text = f"IDOMAIN_L{k + 1}"
            fname = f"idomain.l{k + 1:02d}.bin"
            pth = os.path.join(test.workspace, fname)
            f = open(pth, "wb")
            header = flopy.utils.BinaryHeader.create(
                bintype="HEAD",
                precision="double",
                text=text,
                nrow=nrow,
                ncol=ncol,
                ilay=k + 1,
                pertim=1.0,
                totim=1.0,
                kstp=1,
                kper=1,
            )
            flopy.utils.Util2d.write_bin(
                (nrow, ncol),
                f,
                np.ones((nrow, ncol), dtype=np.int32),
                header_data=header,
            )
            f.close()
            idomain.append(
                {
                    "factor": 1.0,
                    "filename": fname,
                    "data": None,
                    "binary": True,
                    "iprn": 1,
                }
            )
    elif idx == 1:
        fname = "idomain.bin"
        pth = os.path.join(test.workspace, fname)
        f = open(pth, "wb")
        header = flopy.utils.BinaryHeader.create(
            bintype="HEAD",
            precision="double",
            text="IDOMAIN",
            nrow=1,
            ncol=ncol * nrow * nlay,
            ilay=1,
            pertim=1.0,
            totim=1.0,
            kstp=1,
            kper=1,
        )
        flopy.utils.Util2d.write_bin(
            (nrow, ncol),
            f,
            np.ones((nrow * ncol * nlay), dtype=np.int32),
            header_data=header,
        )
        f.close()
        idomain = {
            "factor": 1.0,
            "filename": fname,
            "data": None,
            "binary": True,
            "iprn": 1,
        }

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botarr,
        idomain=idomain,
        filename=f"{name}.dis",
    )

    # initial conditions
    # write initial heads to binary file
    if idx == 0:
        strt = []
        for k in range(nlay):
            text = f"IC_L{k + 1}"
            fname = f"ic.strt_l{k + 1:02d}.bin"
            pth = os.path.join(test.workspace, fname)
            f = open(pth, "wb")
            header = flopy.utils.BinaryHeader.create(
                bintype="HEAD",
                precision="double",
                text=text,
                nrow=nrow,
                ncol=ncol,
                ilay=k + 1,
                pertim=1.0,
                totim=1.0,
                kstp=1,
                kper=1,
            )
            flopy.utils.Util2d.write_bin(
                (nrow, ncol),
                f,
                np.ones((nrow, ncol), dtype=np.float64),
                header_data=header,
            )
            f.close()
            strt.append(
                {
                    "factor": 1.0,
                    "filename": fname,
                    "data": None,
                    "binary": True,
                    "iprn": 1,
                }
            )
    elif idx == 1:
        fname = "ic.strt.bin"
        pth = os.path.join(test.workspace, fname)
        f = open(pth, "wb")
        header = flopy.utils.BinaryHeader.create(
            bintype="HEAD",
            precision="double",
            text="HEAD",
            nrow=1,
            ncol=ncol * nrow * nlay,
            ilay=1,
            pertim=1.0,
            totim=1.0,
            kstp=1,
            kper=1,
        )
        flopy.utils.Util2d.write_bin(
            (nrow, ncol),
            f,
            np.ones((nrow * ncol * nlay), dtype=np.float64),
            header_data=header,
        )
        f.close()
        strt = {
            "factor": 1.0,
            "filename": fname,
            "data": None,
            "binary": True,
            "iprn": 1,
        }

    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    # write icelltype to binary file
    if idx == 0:
        icelltype = []
        for k in range(nlay):
            fname = f"npf.icelltype.l{k + 1}.bin"
            pth = os.path.join(test.workspace, fname)
            f = open(pth, "wb")
            header = flopy.utils.BinaryHeader.create(
                bintype="head",
                text="ICELLTYPE",
                precision="double",
                nrow=nrow,
                ncol=ncol,
                ilay=k + 1,
                pertim=1.0,
                totim=1.0,
                kstp=1,
                kper=1,
            )
            flopy.utils.Util2d.write_bin(
                (nrow, ncol),
                f,
                np.ones((nrow, ncol), dtype=np.int32),
                header_data=header,
            )
            f.close()
            icelltype.append(
                {
                    "factor": 1.0,
                    "filename": fname,
                    "data": None,
                    "binary": True,
                    "iprn": 1,
                }
            )
    elif idx == 1:
        fname = "npf.icelltype.bin"
        pth = os.path.join(test.workspace, fname)
        f = open(pth, "wb")
        header = flopy.utils.BinaryHeader.create(
            bintype="head",
            text="ICELLTYPE",
            precision="double",
            nrow=1,
            ncol=ncol * nrow * nlay,
            ilay=1,
            pertim=1.0,
            totim=1.0,
            kstp=1,
            kper=1,
        )
        flopy.utils.Util2d.write_bin(
            (nrow, ncol),
            f,
            np.ones((nrow * ncol * nlay), dtype=np.int32),
            header_data=header,
        )
        icelltype = {
            "factor": 1.0,
            "filename": fname,
            "data": None,
            "binary": True,
            "iprn": 1,
        }

        f.close()

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=icelltype,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 1.0])
    chdlist0.append([(nlay - 1, nrow - 1, ncol - 1), 0.0])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
    )
    test.run()

"""
Parallel (PETSc) tests for the stress-period-varying IMS linear settings.

Two coupled GWF models are split over two cpus. The period-data file is read
per-rank and applied locally, so no MPI synchronization is required.

Cases:
  - par_ims_period        : default PETSc mode. PERIOD 2 raises INNER_MAXIMUM,
                            switches the preconditioner, and tightens the
                            tolerances; PERIOD 3 lowers INNER_MAXIMUM. Must
                            reproduce the uniform flow field.
  - par_ims_period_badpc  : PRECONDITIONER_LEVELS is rejected when .petscrc
                            selects a native PETSc preconditioner (xfail).
  - par_ims_period_badcnvg: INNER_DVCLOSE is rejected when .petscrc selects a
                            native PETSc convergence check (xfail).
"""

import glob
import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

name_left = "leftmodel"
name_right = "rightmodel"

nper = 3
dis_shape = (5, 5, 5)

# solver data
nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def _inject_period_data(ws, period_text):
    # inject LINEAR_PERIODDATA into the (single) IMS OPTIONS block and write the
    # period-data file; flopy writes "BEGIN options" (mixed case)
    import re

    ims_files = glob.glob(os.path.join(ws, "*.ims"))
    assert len(ims_files) == 1, f"expected one .ims file, found {ims_files}"
    ims = ims_files[0]
    impd = os.path.splitext(os.path.basename(ims))[0] + ".impd"
    with open(ims) as f:
        txt = f.read()
    new = re.sub(
        r"(?i)(BEGIN\s+OPTIONS[^\n]*\n)",
        rf"\1  LINEAR_PERIODDATA FILEIN {impd}\n",
        txt,
        count=1,
    )
    assert new != txt, "could not inject LINEAR_PERIODDATA into the IMS OPTIONS block"
    with open(ims, "w") as f:
        f.write(new)
    with open(os.path.join(ws, impd), "w") as f:
        f.write(period_text)


def get_model(name, ws):
    delr = delc = 100.0
    tops = [0.0, -100.0, -200.0, -300.0, -400.0, -500.0]
    nlay, nrow, ncol = dis_shape
    shift_x = ncol * delr
    h_left, h_right, h_start = 1.0, 10.0, 0.0

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim, time_units="DAYS", nper=nper, perioddata=[(1.0, 1, 1)] * nper
    )
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        relaxation_factor=relax,
        pname="ims",
    )

    for modelname, chd_col, xorigin in (
        (name_left, 0, 0.0),
        (name_right, ncol - 1, shift_x),
    ):
        h = h_left if modelname == name_left else h_right
        gwf = flopy.mf6.ModflowGwf(sim, modelname=modelname, save_flows=True)
        flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            xorigin=xorigin,
            top=tops[0],
            botm=tops[1 : nlay + 1],
        )
        flopy.mf6.ModflowGwfic(gwf, strt=h_start)
        flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=1.0)
        flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data=[
                [(ilay, irow, chd_col), h]
                for irow in range(nrow)
                for ilay in range(nlay)
            ],
        )
        flopy.mf6.ModflowGwfoc(
            gwf,
            head_filerecord=f"{modelname}.hds",
            saverecord=[("HEAD", "LAST")],
        )

    angldegx, cdist = 0.0, delr
    gwfgwf_data = [
        [
            (ilay, irow, ncol - 1),
            (ilay, irow, 0),
            1,
            delr / 2.0,
            delr / 2.0,
            delc,
            angldegx,
            cdist,
        ]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        nexg=len(gwfgwf_data),
        exgmnamea=name_left,
        exgmnameb=name_right,
        exchangedata=gwfgwf_data,
        auxiliary=["ANGLDEGX", "CDIST"],
    )
    return sim


def _write_petscrc(ws, extra_lines):
    with open(os.path.join(ws, ".petscrc"), "w") as f:
        f.write("-ksp_type bcgs\n")
        for line in extra_lines:
            f.write(line + "\n")


def build_default(test):
    sim = get_model("par_ims_period", test.workspace)
    sim.write_simulation()
    _inject_period_data(
        test.workspace,
        "BEGIN PERIOD 2\n"
        "  INNER_MAXIMUM 500\n"
        "  INNER_DVCLOSE 1.0e-8\n"
        "  PRECONDITIONER_LEVELS 5\n"
        "  PRECONDITIONER_DROP_TOLERANCE 1.0e-4\n"
        "END PERIOD\n\n"
        "BEGIN PERIOD 3\n"
        "  INNER_MAXIMUM 100\n"
        "END PERIOD\n",
    )
    return sim, None


def check_default(test):
    fpth = os.path.join(test.workspace, f"{name_left}.hds")
    heads_left = flopy.utils.HeadFile(fpth).get_data().flatten()
    fpth = os.path.join(test.workspace, f"{name_right}.hds")
    heads_right = flopy.utils.HeadFile(fpth).get_data().flatten()
    np.testing.assert_array_almost_equal(heads_left[0:5], [1.0, 2.0, 3.0, 4.0, 5.0])
    np.testing.assert_array_almost_equal(heads_right[0:5], [6.0, 7.0, 8.0, 9.0, 10.0])


def build_reject(name, period_text, petscrc_lines, test):
    sim = get_model(name, test.workspace)
    sim.write_simulation()
    _inject_period_data(test.workspace, period_text)
    _write_petscrc(test.workspace, petscrc_lines)
    return sim, None


@pytest.mark.parallel
def test_default(function_tmpdir, targets):
    test = TestFramework(
        name="par_ims_period",
        workspace=function_tmpdir,
        targets=targets,
        build=build_default,
        check=check_default,
        compare=None,
        parallel=True,
        ncpus=2,
        overwrite=False,  # keep the injected LINEAR_PERIODDATA + period file
    )
    test.run()


reject_cases = [
    (
        "par_ims_period_badpc",
        "BEGIN PERIOD 2\n  PRECONDITIONER_LEVELS 5\nEND PERIOD\n",
        ["-use_petsc_pc", "-sub_pc_type ilu"],
    ),
    (
        "par_ims_period_badcnvg",
        "BEGIN PERIOD 2\n  INNER_DVCLOSE 1.0e-8\nEND PERIOD\n",
        ["-use_petsc_cnvg"],
    ),
]


@pytest.mark.parallel
@pytest.mark.developmode
@pytest.mark.parametrize("name, period_text, petscrc_lines", reject_cases)
def test_reject(name, period_text, petscrc_lines, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_reject(name, period_text, petscrc_lines, t),
        check=None,
        compare=None,
        parallel=True,
        ncpus=2,
        xfail=True,  # the period keyword is rejected under the native PETSc mode
        overwrite=False,
    )
    test.run()

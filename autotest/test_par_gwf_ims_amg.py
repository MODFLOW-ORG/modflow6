"""
Parallel test of the AMG preconditioner (PRECONDITIONER_TYPE AMG) running under
the PETSc solver. A single confined GWF model with a linear head solution is
split into four domains with the flopy Mf6Splitter and run on four ranks. In
parallel the AMG preconditioner is applied per subdomain (block Jacobi), the
same way the MILU0/MILUT preconditioners are, so each rank's local matrix is
solved with its own AMG hierarchy.

Every AMG smoother option is exercised in parallel, one per case:

  - par_amg_v      CG,       default smoother (ILU0 at finest level)
  - par_amg_ilu0a  CG,       SMOOTHER_TYPE ILU0_ALL
  - par_amg_bcgs   BiCGSTAB, default smoother

Each case confirms that (1) the parallel AMG solve reproduces the exact linear
solution, and (2) the requested AMG sub-preconditioner options are actually
active (the PETSc solver reports "IMS AMG" plus the smoother type), so a silent
fall-back to ILU or to the default options cannot pass unnoticed.

Adaptive omega (RELAXATION_FACTOR 0) is a variable preconditioner and is not
allowed under a Krylov accelerator, so the AMG cases use RELAXATION_FACTOR 1.
"""

from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "par_amg_v",  # default smoother (ILU0 at finest), CG
    "par_amg_ilu0a",  # SMOOTHER_TYPE ILU0_ALL, CG
    "par_amg_bcgs",  # default smoother, BiCGSTAB
]
# per case: linear_acceleration, smoother_type, expected smoother label
amg_opts = [
    ("CG", None, "ILU0"),
    ("CG", "ILU0_ALL", "ILU0_ALL"),
    ("BICGSTAB", None, "ILU0"),
]

gwf_name = "gwf"

# solver data
nouter, ninner = 100, 300
hclose, rclose = 1.0e-9, 1.0e-3
amg_levels = 4

# boundary heads (linear solution between them)
h_left = 1.0
h_right = 10.0

# model discretization
nlay, nrow, ncol = 1, 20, 20
delr = delc = 50.0


def _split4(sim, test, smoother_type=None):
    """Split a single-model simulation into four domains, one per rank, and
    re-assert the AMG options on the split IMS (defensive, in case the splitter
    regenerates the solution with defaults)."""
    mfsplit = flopy.mf6.utils.Mf6Splitter(sim)
    split_array = np.zeros((nrow, ncol), dtype=int)
    for irow in range(nrow):
        for icol in range(ncol):
            idomain = 0
            if irow > nrow / 2:
                idomain += 2
            if icol > ncol / 2:
                idomain += 1
            split_array[irow, icol] = idomain
    split_sim = mfsplit.split_model(split_array)
    split_sim.set_sim_path(test.workspace)

    for pkg in split_sim.sim_package_list:
        if pkg.package_type.lower() == "ims":
            pkg.preconditioner_type = "AMG"
            pkg.relaxation_factor = 1.0
            pkg.preconditioner_levels = amg_levels
            pkg.smoother_type = smoother_type

    partitions = [(m, i) for i, m in enumerate(split_sim.model_names)]
    flopy.mf6.ModflowUtlhpc(split_sim, partitions=partitions)
    return split_sim


def get_model(idx, test):
    name = cases[idx]
    accel, smoother_type, _ = amg_opts[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=test.workspace
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])

    # AMG preconditioner; RELAXATION_FACTOR must be 1 in parallel (Krylov needs a
    # fixed preconditioner, so adaptive omega is not supported)
    flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=accel,
        preconditioner_type="AMG",
        relaxation_factor=1.0,
        preconditioner_levels=amg_levels,
        smoother_type=smoother_type,
        pname="ims",
    )

    left_chd = [
        [(ilay, irow, 0), h_left] for irow in range(nrow) for ilay in range(nlay)
    ]
    right_chd = [
        [(ilay, irow, ncol - 1), h_right]
        for irow in range(nrow)
        for ilay in range(nlay)
    ]
    chd_spd = {0: left_chd + right_chd}

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwf_name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=[-10.0],
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{gwf_name}.hds",
        saverecord=[("HEAD", "LAST")],
    )

    return _split4(sim, test, smoother_type=smoother_type)


def build_models(idx, test):
    return get_model(idx, test), None


def check_output(idx, test):
    _, _, exp_smoother = amg_opts[idx]

    # exact linear head between the two constant-head boundaries
    def exact(x):
        return h_left + (h_right - h_left) * (x - 0.5 * delr) / ((ncol - 1) * delr)

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)

    # every sub-model must reproduce the exact solution
    for model_name in sim.model_names:
        hds = flopy.utils.HeadFile(
            Path(test.workspace) / f"{model_name}.hds"
        ).get_data()
        grb = flopy.mf6.utils.MfGrdFile(Path(test.workspace) / f"{model_name}.dis.grb")
        xyc = grb.modelgrid.xycenters
        xoff = grb.modelgrid.xoffset
        for ilay in range(grb.modelgrid.nlay):
            for irow in range(grb.modelgrid.nrow):
                for icol in range(grb.modelgrid.ncol):
                    xc = xyc[0][icol] + xoff
                    assert hds[ilay, irow, icol] == pytest.approx(
                        exact(xc), abs=1e-4
                    ), f"model {model_name} cell ({ilay},{irow},{icol}) head mismatch"

    # confirm AMG (and the requested smoother) was actually active; the PETSc
    # solver prints "IMS AMG" and the AMG smoother type
    listings = sorted(Path(test.workspace).glob("*.lst"))
    text = "\n".join(open(p).read() for p in listings)
    assert "IMS AMG" in text, "PETSc solver did not report the AMG sub-preconditioner"

    smoother_lines = [ln for ln in text.splitlines() if "AMG smoother type" in ln]
    assert smoother_lines and all(exp_smoother in ln for ln in smoother_lines), (
        f"expected AMG smoother {exp_smoother} in the listing"
    )


@pytest.mark.parallel
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        compare=None,
        parallel=True,
        ncpus=4,
    )
    test.run()


@pytest.mark.parallel
def test_amg_relax0_error(function_tmpdir, targets):
    """AMG in parallel with RELAXATION_FACTOR 0 (adaptive omega) must terminate
    with an error: a variable preconditioner is incompatible with the PETSc
    Krylov accelerator."""

    def build(test):
        split = get_model(0, test)  # CG + AMG, RELAXATION_FACTOR 1
        for pkg in split.sim_package_list:
            if pkg.package_type.lower() == "ims":
                pkg.relaxation_factor = 0.0
        split.write_simulation(silent=True)
        return split

    test = TestFramework(
        name="par_amg_r0_err",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=None,
        compare=None,
        overwrite=False,
        parallel=True,
        ncpus=4,
        xfail=True,
    )
    test.run()


# Newton (unconfined) coverage: areal 1-D unconfined flow driven by west/east
# constant heads. The Newton Jacobian is nonsymmetric (BiCGSTAB) and its values
# change every outer iteration, which exercises the AMG value-update path in
# parallel. Correctness is checked against a serial ILU run of the same split
# model (identical discrete system, independent solver) rather than an analytic
# reference, so the comparison is tight and free of discretization bias.
h_west_nwt = 9.0
h_east_nwt = 2.0


def get_newton_model(test, sim_ws, amg, hpc):
    sim = flopy.mf6.MFSimulation(
        sim_name="par_amg_nwt", version="mf6", exe_name="mf6", sim_ws=str(sim_ws)
    )
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)])
    ims_kw = dict(
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",  # Newton Jacobian is nonsymmetric
        pname="ims",
    )
    if amg:
        ims_kw.update(
            preconditioner_type="AMG",
            relaxation_factor=1.0,
            preconditioner_levels=amg_levels,
        )
    flopy.mf6.ModflowIms(sim, **ims_kw)
    left_chd = [[(0, irow, 0), h_west_nwt] for irow in range(nrow)]
    right_chd = [[(0, irow, ncol - 1), h_east_nwt] for irow in range(nrow)]
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwf_name, newtonoptions="NEWTON", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=nrow, ncol=ncol, delr=delr, delc=delc, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfic(gwf, strt=h_west_nwt)
    flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=10.0)  # unconfined
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: left_chd + right_chd})
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{gwf_name}.hds", saverecord=[("HEAD", "LAST")]
    )

    mfsplit = flopy.mf6.utils.Mf6Splitter(sim)
    split_array = np.zeros((nrow, ncol), dtype=int)
    for irow in range(nrow):
        for icol in range(ncol):
            split_array[irow, icol] = (2 if irow > nrow / 2 else 0) + (
                1 if icol > ncol / 2 else 0
            )
    split_sim = mfsplit.split_model(split_array)
    split_sim.set_sim_path(str(sim_ws))
    if amg:
        for pkg in split_sim.sim_package_list:
            if pkg.package_type.lower() == "ims":
                pkg.preconditioner_type = "AMG"
                pkg.relaxation_factor = 1.0
                pkg.preconditioner_levels = amg_levels
    if hpc:
        partitions = [(m, i) for i, m in enumerate(split_sim.model_names)]
        flopy.mf6.ModflowUtlhpc(split_sim, partitions=partitions)
    return split_sim


def check_newton(test):
    ser_ws = Path(test.workspace) / "serial"
    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    for model_name in sim.model_names:
        h_par = flopy.utils.HeadFile(
            Path(test.workspace) / f"{model_name}.hds"
        ).get_data()
        h_ser = flopy.utils.HeadFile(ser_ws / f"{model_name}.hds").get_data()
        maxdiff = float(np.nanmax(np.abs(h_par - h_ser)))
        assert maxdiff < 1e-4, (
            f"parallel AMG vs serial ILU head mismatch for {model_name}: {maxdiff}"
        )

    listings = sorted(Path(test.workspace).glob("*.lst"))
    assert any("IMS AMG" in open(p).read() for p in listings), "AMG not active"


@pytest.mark.parallel
def test_amg_newton(function_tmpdir, targets):
    def build(test):
        parallel = get_newton_model(test, test.workspace, amg=True, hpc=True)
        serial = get_newton_model(
            test, Path(test.workspace) / "serial", amg=False, hpc=False
        )
        return parallel, serial

    test = TestFramework(
        name="par_amg_nwt",
        workspace=function_tmpdir,
        targets=targets,
        build=build,
        check=check_newton,
        compare=None,
        parallel=True,
        ncpus=[4, 1],
    )
    test.run()

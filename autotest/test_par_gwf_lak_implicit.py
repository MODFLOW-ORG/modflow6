"""
Parallel integration test for the LAK IMPLICIT option, based on Lake Package
problem 2 (ex-gwf-lak-p02; Merritt and Konikow, 2000).

A 5x27x17 model with two lakes -- a northern lake and a southern lake, each a
5x5 (layer 0) / 3x3 (layer 1) bowl -- separated in the row (y) direction with a
gap between them, plus SFR streamflow routing that runs north to south through
the domain. The two lakes start in a single LAK package.

The flopy model splitter divides the domain into two halves at y ~ 10000 ft,
which falls in the gap between the lakes, so each lake ends up entirely in one
sub-model (and therefore in its own LAK package). The splitter sets up the
inter-model SFR mover for the stream reaches that cross the split and the
GWF-GWF exchange between the two sub-models. The IMPLICIT option is then enabled
on the northern lake.

The split simulation is run two ways: once on two processes (the main run) and
once serially (the "mf6" comparison run). The test framework compares their
heads automatically (compare="mf6"), and the check function additionally
compares their cell budgets. A lake is a single dependent variable that lives
entirely on one process, so it is never split across domains -- this test confirms
that the implicit lake formulation (which makes the coefficient matrix asymmetric
and requires BiCGSTAB) is assembled and solved correctly by the parallel (PETSc)
solver, giving the same result as the serial run.

The example's lake-to-stream movers (SFR <-> LAK both ways for each lake) are
added back manually after the split: the flopy splitter cannot remap a
model-level MVR across a split, but these movers are intra-domain here -- each
lake and its stream reaches stay in the same sub-model -- so each is added as a
model-level MVR with the reaches renumbered to the sub-model's local indices.
"""

import flopy
import numpy as np
import pytest
from flopy.mf6.utils import Mf6Splitter, get_lak_connections
from flopy.utils.compare import compare_cell_budget
from framework import TestFramework

# model parameters (from ex-gwf-lak-p02)
name = "lakp02"
nlay, nrow, ncol = 5, 27, 17
top = 200.0
botm = [102.0, 97.0, 87.0, 77.0, 67.0]
strt = 115.0
k11 = k33 = 30.0
ss, sy = 3e-4, 0.2
H1, H2 = 160.0, 140.0
recharge = 0.0116
etvrate, etvdepth = 0.0141, 15.0
lak_strt, lak_etrate, lak_bedleak = 130.0, 0.0103, 0.1
# reduced from the published 200 steps to keep the parallel test fast
tdis_ds = ((1500.0, 100, 1.02),)
y_split = 10000.0

delr = np.array(
    [250.0, 1000, 1000, 1000, 1000, 1000, 500, 500, 500, 500, 500,
     1000, 1000, 1000, 1000, 1000, 250.0]
)  # fmt: skip
delc = np.array(
    [250.0, 1000, 1000, 1000, 1000, 1000, 500, 500, 500, 500, 500,
     1000, 1000, 1000, 1000, 1000, 500, 500, 500, 500, 500,
     1000, 1000, 1000, 1000, 1000, 250.0]
)  # fmt: skip
shape3d = (nlay, nrow, ncol)

lak_outlets = [
    [0, 0, -1, "manning", 114.85, 5.0, 0.05, 8.206324419006205e-4],
    [1, 1, -1, "manning", 109.4286, 5.0, 0.05, 9.458197164349258e-4],
]
lak_spd = [
    [0, "rainfall", recharge],
    [0, "evaporation", lak_etrate],
    [1, "rainfall", recharge],
    [1, "evaporation", lak_etrate],
]

# SFR routes north (reach 0) to south (reach 21), crossing the y=10000 split
sfr_pakdata = [
    [0, 0, 1, 4, 1000, 5, 0.001103448, 123.94827, 0.5, 0.5, 0.050000001, 1, 1, 0],
    [1, 0, 2, 4, 1000, 5, 0.001103448, 122.84483, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [2, 0, 3, 4, 1000, 5, 0.001103448, 121.74138, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [3, 0, 3, 5, 1000, 5, 0.001103448, 120.63793, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [4, 0, 3, 6, 500, 5, 0.001103448, 119.81035, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [5, 0, 3, 7, 750, 5, 0.001103448, 119.12069, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [6, 0, 4, 7, 1000, 5, 0.001103448, 118.15517, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [7, 0, 5, 7, 1000, 5, 0.001103448, 117.05173, 0.5, 0.5, 0.050000001, 1, 1, 0],
    [8, 0, 11, 8, 1000, 5, 0.000820632, 114.43968, 0.5, 0.5, 0.050000001, 1, 1, 0],
    [9, 0, 12, 8, 1000, 5, 0.000820632, 113.61905, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [10, 0, 13, 9, 559, 5, 0.000820632, 112.97937, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [11, 0, 13, 9, 559, 5, 0.000820632, 112.52063, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [12, 0, 14, 9, 1000, 5, 0.000820632, 111.88095, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [13, 0, 15, 9, 1000, 5, 0.000820632, 111.06032, 0.5, 0.5, 0.050000001, 1, 1, 0],
    [14, 0, 21, 9, 1000, 5, 0.00094582, 108.95569, 0.5, 0.5, 0.050000001, 1, 1, 0],
    [15, 0, 22, 9, 750, 5, 0.00094582, 108.1281, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [16, 0, 22, 10, 500, 5, 0.00094582, 107.53696, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [17, 0, 22, 11, 1000, 5, 0.00094582, 106.82759, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [18, 0, 22, 12, 1000, 5, 0.00094582, 105.88177, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [19, 0, 22, 13, 1000, 5, 0.00094582, 104.93595, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [20, 0, 22, 14, 1000, 5, 0.00094582, 103.99014, 0.5, 0.5, 0.050000001, 2, 1, 0],
    [21, 0, 22, 15, 1000, 5, 0.00094582, 103.04431, 0.5, 0.5, 0.050000001, 1, 1, 0],
]
sfr_conn = [
    [0, -1], [1, 0, -2], [2, 1, -3], [3, 2, -4], [4, 3, -5], [5, 4, -6],
    [6, 5, -7], [7, 6], [8, -9], [9, 8, -10], [10, 9, -11], [11, 10, -12],
    [12, 11, -13], [13, 12], [14, -15], [15, 14, -16], [16, 15, -17],
    [17, 16, -18], [18, 17, -19], [19, 18, -20], [20, 19, -21], [21, 20],
]  # fmt: skip
sfr_spd = [[0, "inflow", 691200.0]]


def _lake_map():
    # northern lake (index 0) and southern lake (index 1), each a 5x5 footprint
    # in layer 0 with a 3x3 inner cell in layer 1, separated by a gap of rows
    lm = -np.ones(shape3d, dtype=int)
    lm[0, 6:11, 6:11] = 0
    lm[1, 7:10, 7:10] = 0
    lm[0, 16:21, 6:11] = 1
    lm[1, 17:20, 7:10] = 1
    return lm


def _et_surface(lake_map):
    xlen = delr.sum() - 0.5 * (delr[0] + delr[-1])
    x = 0.0
    s1d = H1 * np.ones(ncol)
    for i in range(1, ncol):
        x += 0.5 * (delr[i - 1] + delr[i])
        s1d[i] = H1 + (H2 - H1) * (x / xlen)
    surf = np.tile(s1d, (nrow, 1))
    surf[lake_map[0] > -1] = botm[0] - 2
    surf[lake_map[1] > -1] = botm[1] - 2
    return surf


def _build_base(ws, exe):
    lake_map = _lake_map()
    chd_spd = []
    for k in range(nlay):
        chd_spd += [[k, i, 0, H1] for i in range(nrow)]
        chd_spd += [[k, i, ncol - 1, H2] for i in range(nrow)]

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=tdis_ds, time_units="days")
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        linear_acceleration="bicgstab",
        outer_maximum=500,
        outer_dvclose=1e-9,
        inner_maximum=100,
        inner_dvclose=1e-9,
        rcloserecord="1e-6 strict",
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="newton", save_flows=True
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units="feet",
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        idomain=np.ones(shape3d, dtype=int),
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k11, k33=k33)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=1, sy=sy, ss=ss)
    flopy.mf6.ModflowGwfic(gwf, strt=strt)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd_spd)
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    flopy.mf6.ModflowGwfevta(
        gwf, surface=_et_surface(lake_map), rate=etvrate, depth=etvdepth
    )
    idomain_wlakes, pakdata_dict, lak_conn = get_lak_connections(
        gwf.modelgrid, lake_map, bedleak=lak_bedleak
    )
    lak_packagedata = [[k, lak_strt, pakdata_dict[k]] for k in pakdata_dict]
    flopy.mf6.ModflowGwflak(
        gwf,
        pname="LAK-1",
        time_conversion=86400.0,
        length_conversion=3.28081,
        mover=True,
        print_stage=True,
        nlakes=2,
        noutlets=len(lak_outlets),
        packagedata=lak_packagedata,
        connectiondata=lak_conn,
        outlets=lak_outlets,
        perioddata=lak_spd,
    )
    gwf.dis.idomain = idomain_wlakes
    flopy.mf6.ModflowGwfsfr(
        gwf,
        pname="SFR-1",
        time_conversion=86400.0,
        length_conversion=3.28081,
        mover=True,
        print_stage=True,
        nreaches=len(sfr_pakdata),
        packagedata=sfr_pakdata,
        connectiondata=sfr_conn,
        perioddata=sfr_spd,
    )
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )
    return sim


def _lak_package(model):
    # return the model's LAK package by type rather than by the splitter's
    # generated package name (which is an implementation detail of Mf6Splitter)
    for pkg in model.packagelist:
        if pkg.package_type.lower() == "lak":
            return pkg
    raise KeyError(f"no LAK package found in model {model.name}")


def _split_with_movers(ws, exe, hpc):
    # build the base model, split it at y ~ 10000 ft, add the intra-domain
    # lake<->stream movers the splitter cannot remap, enable IMPLICIT on the
    # northern lake, and optionally add the HPC partitioning
    base = _build_base(str(ws), exe)
    gwf = base.get_model(name)
    split_array = np.where(gwf.modelgrid.ycellcenters >= y_split, 0, 1).astype(int)

    # each lake must end up entirely within one sub-domain
    lake_map = _lake_map()
    for li in (0, 1):
        doms = set(split_array[lake_map[0] == li].tolist())
        doms |= set(split_array[lake_map[1] == li].tolist())
        assert len(doms) == 1, f"lake {li} straddles the split: {doms}"

    split = Mf6Splitter(base).split_model(split_array)
    north = f"{name}_0"  # sub-model for the northern lake (domain 0, y >= split)
    south = f"{name}_1"  # sub-model for the southern lake
    assert north in split.model_names and south in split.model_names

    # The example's lake<->stream movers are intra-domain after the split -- each
    # lake and its stream reaches stay in the same sub-model -- so add each as a
    # model-level MVR. Reaches keep their original 0-based numbering in the
    # northern sub-model (reaches 7 and 8); in the southern sub-model the original
    # reaches 13 and 14 become local reaches 1 and 2.
    flopy.mf6.ModflowGwfmvr(
        split.get_model(north),
        maxmvr=2,
        maxpackages=2,
        packages=[("SFR-1",), ("LAK-1",)],
        perioddata=[
            ["SFR-1", 7, "LAK-1", 0, "FACTOR", 1.0],
            ["LAK-1", 0, "SFR-1", 8, "FACTOR", 1.0],
        ],
    )
    flopy.mf6.ModflowGwfmvr(
        split.get_model(south),
        maxmvr=2,
        maxpackages=2,
        packages=[("SFR-1",), ("LAK-1",)],
        perioddata=[
            ["SFR-1", 1, "LAK-1", 0, "FACTOR", 1.0],
            ["LAK-1", 0, "SFR-1", 2, "FACTOR", 0.5],
        ],
    )

    # enable the IMPLICIT formulation on the northern lake
    _lak_package(split.get_model(north)).implicit = True

    if hpc:
        partitions = [(m, i) for i, m in enumerate(split.model_names)]
        flopy.mf6.ModflowUtlhpc(split, partitions=partitions)

    split.set_sim_path(str(ws))
    return split


def build_models(test):
    # Build the split simulation twice: the first (in the workspace) is run in
    # parallel and the second (in the "mf6" subdirectory) serially as the
    # comparison; the framework writes both and compares their heads and budgets
    # (compare="mf6").
    exe = test.targets["mf6"]
    sim = _split_with_movers(test.workspace, exe, hpc=True)
    cmp = _split_with_movers(test.workspace / "mf6", exe, hpc=False)
    return sim, cmp


def check_output(test):
    # confirm the northern lake actually used the implicit formulation in the
    # parallel run
    nlst = (test.workspace / f"{name}_0.lst").read_text()
    assert "IMPLICIT FORMULATION" in nlst, "northern lake did not use IMPLICIT"

    # the framework comparison (compare="mf6") checks the heads of the parallel
    # (workspace) and serial ("mf6" subdirectory) runs; also compare their cell
    # budgets here
    for model in (f"{name}_0", f"{name}_1"):
        assert compare_cell_budget(
            test.workspace / f"{model}.cbc",
            test.workspace / "mf6" / f"{model}.cbc",
            outfile=test.workspace / f"{model}.cbc.cmp",
            rclose=1e-3,
        ), f"serial vs parallel budget mismatch for {model}"


@pytest.mark.parallel
def test_mf6model(function_tmpdir, targets):
    test = TestFramework(
        name="par_lak_implicit",
        workspace=function_tmpdir,
        targets=targets,
        build=build_models,
        check=check_output,
        compare="mf6",
        parallel=True,
        # build_models returns two simulations: index 0 is the split model run in
        # parallel on 2 ranks, index 1 is the "mf6" comparison run serially on 1
        ncpus=[2, 1],
    )
    test.run()

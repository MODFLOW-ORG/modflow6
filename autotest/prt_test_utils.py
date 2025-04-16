import os

import flopy
import matplotlib as mpl
import numpy as np

DEFAULT_EXIT_SOLVE_TOL = 1e-5


def all_equal(series, val):
    a = series.to_numpy()
    return a[0] == val and (a[0] == a).all()


class HorizontalCase:
    nlay = 1
    nrow = 1
    ncol = 10
    top = 1
    botm = [0.0]
    nper = 1
    perlen = 1.0
    nstp = 1
    tsmult = 1.0
    porosity = 0.1
    releasepts_mp7 = [
        # node number, localx, localy, localz
        (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
        for i in range(9)
    ]
    releasepts_prt = [
        # particle index, k, i, j, x, y, z
        [i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5]
        for i in range(9)
    ]

    @staticmethod
    def get_gwf_sim(name, ws, mf6) -> flopy.mf6.MFSimulation:
        """
        Simple GWF simulation on a simple horizontal line grid.
        """

        # create simulation
        sim = flopy.mf6.MFSimulation(
            sim_name=name,
            exe_name=mf6,
            version="mf6",
            sim_ws=ws,
        )

        # create tdis package
        flopy.mf6.modflow.mftdis.ModflowTdis(
            sim,
            pname="tdis",
            time_units="DAYS",
            nper=HorizontalCase.nper,
            perioddata=[
                (
                    HorizontalCase.perlen,
                    HorizontalCase.nstp,
                    HorizontalCase.tsmult,
                )
            ],
        )

        # create gwf model
        gwfname = f"{name}_gwf"
        gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

        # create gwf discretization
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            gwf,
            pname="dis",
            nlay=HorizontalCase.nlay,
            nrow=HorizontalCase.nrow,
            ncol=HorizontalCase.ncol,
        )

        # create gwf initial conditions package
        flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

        # create gwf node property flow package
        flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            gwf,
            pname="npf",
            save_saturation=True,
            save_specific_discharge=True,
        )

        # create gwf chd package
        spd = {
            0: [[(0, 0, 0), 1.0, 1.0], [(0, 0, 9), 0.0, 0.0]],
            1: [[(0, 0, 0), 0.0, 0.0], [(0, 0, 9), 1.0, 2.0]],
        }
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            pname="CHD-1",
            stress_period_data=spd,
            auxiliary=["concentration"],
        )

        # create gwf output control package
        # output file names
        gwf_budget_file = f"{gwfname}.bud"
        gwf_head_file = f"{gwfname}.hds"
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=gwf_budget_file,
            head_filerecord=gwf_head_file,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # create iterative model solution for gwf model
        ims = flopy.mf6.ModflowIms(sim)

        return sim


class FlopyReadmeCase:
    nlay = 1
    nrow = 10
    ncol = 10
    top = 1.0
    botm = [0.0]
    nper = 1
    perlen = 1.0
    nstp = 1
    tsmult = 1.0
    porosity = 0.1
    releasepts_mp7 = [
        # node number, localx, localy, localz
        (0, float(f"0.{i + 1}"), float(f"0.{i + 1}"), 0.5)
        for i in range(9)
    ]
    releasepts_prt = [
        # particle index, k, i, j, x, y, z
        [i, 0, 0, 0, float(f"0.{i + 1}"), float(f"9.{i + 1}"), 0.5]
        for i in range(9)
    ]

    @staticmethod
    def get_gwf_sim(name, ws, mf6) -> flopy.mf6.MFSimulation:
        """
        Simple GWF simulation for use/modification by PRT tests
        """

        # create simulation
        sim = flopy.mf6.MFSimulation(
            sim_name=name,
            exe_name=mf6,
            version="mf6",
            sim_ws=ws,
        )

        # create tdis package
        flopy.mf6.modflow.mftdis.ModflowTdis(
            sim,
            pname="tdis",
            time_units="DAYS",
            nper=FlopyReadmeCase.nper,
            perioddata=[
                (
                    FlopyReadmeCase.perlen,
                    FlopyReadmeCase.nstp,
                    FlopyReadmeCase.tsmult,
                )
            ],
        )

        # create gwf model
        gwfname = f"{name}_gwf"
        gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

        # create gwf discretization
        flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(
            gwf,
            pname="dis",
            nlay=FlopyReadmeCase.nlay,
            nrow=FlopyReadmeCase.nrow,
            ncol=FlopyReadmeCase.ncol,
            top=FlopyReadmeCase.top,
            botm=FlopyReadmeCase.botm,
        )

        # create gwf initial conditions package
        flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, pname="ic")

        # create gwf node property flow package
        flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(
            gwf,
            pname="npf",
            save_saturation=True,
            save_specific_discharge=True,
        )

        # create gwf chd package
        spd = {
            0: [[(0, 0, 0), 1.0, 1.0], [(0, 9, 9), 0.0, 0.0]],
            1: [[(0, 0, 0), 0.0, 0.0], [(0, 9, 9), 1.0, 2.0]],
        }
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            pname="CHD-1",
            stress_period_data=spd,
            auxiliary=["concentration"],
        )

        # create gwf output control package
        # output file names
        gwf_budget_file = f"{gwfname}.bud"
        gwf_head_file = f"{gwfname}.hds"
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=gwf_budget_file,
            head_filerecord=gwf_head_file,
            saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        )

        # create iterative model solution for gwf model
        ims = flopy.mf6.ModflowIms(sim)

        return sim


def check_track_data(
    track_bin: os.PathLike,
    track_hdr: os.PathLike,
    track_csv: os.PathLike,
):
    """Check that binary and CSV track files are equal."""

    # get dtype from ascii header file
    dt = get_track_dtype(track_hdr)

    # read output files
    data_bin = np.fromfile(track_bin, dtype=dt)
    data_csv = np.genfromtxt(track_csv, dtype=dt, delimiter=",", names=True)
    if len(data_csv.shape) == 0:
        # https://stackoverflow.com/a/24943993/6514033
        data_csv = np.array([data_csv])

    # check shape
    assert data_bin.shape == data_csv.shape, (
        "Binary and CSV track data shapes do not match: "
        f"{data_bin.shape} != {data_csv.shape}"
    )

    # check each column separately to avoid TypeError:
    # The DType <class 'numpy._FloatAbstractDType'> could not be promoted by <class 'numpy.dtype[void]'>  # noqa
    for k in data_bin.dtype.names:
        if k == "name":
            continue
        assert np.allclose(data_bin[k], data_csv[k], equal_nan=True)

    # make sure columns have values in the expected range
    assert all(data_bin["iprp"] >= 1)
    assert all(data_bin["irpt"] >= 1)
    assert all(data_bin["kper"] >= 1)
    assert all(data_bin["kstp"] >= 1)
    assert all(data_bin["ilay"] >= 1)
    assert all(data_bin["icell"] >= 1)
    assert all(data_bin["istatus"] >= 0)
    assert all(data_bin["ireason"] >= 0)


def check_budget_data(lst: os.PathLike, perlen=1, nper=1, nstp=1):
    # load PRT model's list file
    mflist = flopy.utils.mflistfile.ListBudget(
        lst, budgetkey="MASS BUDGET FOR ENTIRE MODEL"
    )
    names = mflist.get_record_names()
    entries = mflist.entries

    # check timesteps
    inc = mflist.get_incremental()
    v = inc["totim"][-1]
    exp = float(perlen * nper)
    assert v == exp, f"Last time should be {exp}, found {v}"

    # entries should be a subset of names
    assert all(e in names for e in entries)

    expected_entries = [
        "PRP_IN",
        "PRP_OUT",
    ]
    assert all(en in names for en in expected_entries)


def get_model_name(name, mdl):
    return f"{name}_{mdl}"


def get_track_dtype(path: os.PathLike):
    """Read a numpy dtype describing particle track
    data format from the ascii track header file."""

    hdr_lns = open(path).readlines()
    hdr_lns_spl = [[ll.strip() for ll in l.split(",")] for l in hdr_lns]
    return np.dtype(list(zip(hdr_lns_spl[0], hdr_lns_spl[1])))


def get_ireason_code(output_event):
    """
    Map output event to PRT ireason code specifying
    the reason a particle track datum was recorded.
    """

    return (
        0
        if output_event == "RELEASE"
        else (
            1
            if output_event == "TRANSIT"
            else (
                2
                if output_event == "TIMESTEP"
                else (
                    3
                    if output_event == "TERMINATE"
                    else 4
                    if output_event == "WEAKSINK"
                    else -1
                )
            )
        )
    )


def get_partdata(grid, rpts):
    """
    Make a flopy.modpath.ParticleData from the given grid and release points.
    """

    if grid.grid_type == "structured":
        return flopy.modpath.ParticleData(
            partlocs=[grid.get_lrc(p[0])[0] for p in rpts],
            structured=True,
            localx=[p[1] for p in rpts],
            localy=[p[2] for p in rpts],
            localz=[p[3] for p in rpts],
            timeoffset=0,
            drape=0,
        )
    else:
        return flopy.modpath.ParticleData(
            partlocs=[p[0] for p in rpts],
            structured=False,
            localx=[p[1] for p in rpts],
            localy=[p[2] for p in rpts],
            localz=[p[3] for p in rpts],
            timeoffset=0,
            drape=0,
        )


def has_default_boundnames(data):
    name = [int(n.partition("0")[2]) for n in data["name"].to_numpy()]
    irpt = data["irpt"].to_numpy()
    return np.array_equal(name, irpt)


def plot_nodes_and_vertices(
    gwf, mg, ibd, ncpl, ax, xmin=None, xmax=None, ymin=None, ymax=None
):
    """
    Plot cell nodes and vertices (and IDs) on a zoomed inset
    """

    ax.set_aspect("equal")
    xlim = False
    ylim = False
    if xmin is not None and xmax is not None:
        ax.set_xlim([xmin, xmax])
        xlim = True
    if ymin is not None and ymax is not None:
        ax.set_ylim([ymin, ymax])
        ylim = True

    # create map view plot
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    v = pmv.plot_grid(lw=0.5, edgecolor="black")
    t = ax.set_title("Node and vertex indices (one-based)\n", fontsize=14)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])

    # plot vertices
    verts = mg.verts
    ax.plot(verts[:, 0], verts[:, 1], "bo", alpha=0.5)
    for i in range(ncpl):
        x, y = verts[i, 0], verts[i, 1]
        ax.annotate(str(i + 1), verts[i, :], color="b", alpha=0.5)

    # plot nodes
    xc, yc = mg.get_xcellcenters_for_layer(0), mg.get_ycellcenters_for_layer(0)
    if mg.grid_type == "structured":
        ids = mg.get_lrc(range(mg.ncpl))
        for _, i, j in ids:
            x, y = xc[i], yc[i]
            ax.plot(x, y, "o", color="grey", alpha=0.5)
            ax.annotate(f"{i + 1}, {j + 1}", (x, y), color="grey", alpha=0.5)
    else:
        for i in range(mg.ncpl):
            x, y = xc[i], yc[i]
            ax.plot(x, y, "o", color="grey", alpha=0.5)
            ax.annotate(str(i + 1), (x, y), color="grey", alpha=0.5)

    # create legend
    ax.legend(
        handles=[
            mpl.patches.Patch(color="blue", label="vertex"),
            mpl.patches.Patch(color="grey", label="node"),
        ],
        loc="upper left",
    )

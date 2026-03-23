"""
Test local z coordinate conversion for PRT PRP package.

This test verifies two fixes in the particle local-to-model z coordinate conversion
when using the LOCALZ option in the PRP package:

Fix #1: Constrain effective top to cell vertical extent
    - For convertible cells, whose saturated thickness depends on head, set effective
      cell top used for local z conversion no higher than the geometric cell top and
      no lower than the cell bottom.

Fix #2: Use correct top based on cell type
    - Confined cells (icelltype==0): Use geometric top for local z conversion
    - Convertible cells (icelltype!=0): Use head (subject to constraint in fix #1) as
      effective top for local z conversion

The test includes 6 cases, each with a simple single-layer model:
1. Convertible cell with head within bounds → uses head as top
2. Confined cell with head within bounds → uses geometric top
3. Convertible cell with head < bottom (dry cell) → clamps to bottom
4. Confined cell with head < bottom → uses geometric top (ignores head)
5. Convertible cell with head > top → clamps to top
6. Confined cell with head > top → clamps to top

Each case has both FMI (separate sims) and exchange (same sim) variants.

Local z coordinate conversion formula:
- Confined: z_model = bot + z_local * (top_geometric - bot)
- Convertible: z_model = bot + z_local * (max(head, bot) - bot)
"""

from pathlib import Path

import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from prt_test_utils import get_model_name

simname = "prtz"

# Six test cases, each with FMI and exchange variants
cases = [
    f"{simname}cvt",  # convertible, bottom < head < top
    f"{simname}cvtexg",
    f"{simname}cnf",  # confined, bottom < head < top
    f"{simname}cnfexg",
    f"{simname}cvtdr",  # convertible, head < bottom
    f"{simname}cvtexgdr",
    f"{simname}cnflow",  # confined, head < bottom
    f"{simname}cnfexglo",
    f"{simname}cvthi",  # convertible, head > top
    f"{simname}cvtexghi",
    f"{simname}cnfhi",  # confined, head > top
    f"{simname}cnfexghi",
]

# Model parameters (same for all cases)
nlay, nrow, ncol = 1, 1, 1
delr = delc = 10.0
nper = 1
perlen = 1.0
nstp = 1
tsmult = 1.0
porosity = 0.1

# Case-specific parameters
case_params = {
    "cvt": {  # Case 1: Convertible, head within bounds
        "top": 10.0,
        "botm": [5.0],
        "strt": [[7.5]],
        "icelltype": [1],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            0: 5.0,  # bot + 0.0 * (head - bot) = 5.0 + 0.0 * 2.5 = 5.0
            1: 6.25,  # bot + 0.5 * (head - bot) = 5.0 + 0.5 * 2.5 = 6.25
            2: 7.5,  # bot + 1.0 * (head - bot) = 5.0 + 1.0 * 2.5 = 7.5 (head)
        },
        "description": "Convertible cell, head within bounds → head as top",
    },
    "cnf": {  # Case 2: Confined, head within bounds
        "top": 10.0,
        "botm": [5.0],
        "strt": [[7.5]],
        "icelltype": [0],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            0: 5.0,  # bot + 0.0 * (top - bot) = 5.0 + 0.0 * 5.0 = 5.0
            1: 7.5,  # bot + 0.5 * (top - bot) = 5.0 + 0.5 * 5.0 = 7.5
            2: 10.0,  # bot + 1.0 * (top - bot) = 5.0 + 1.0 * 5.0 = 10.0 (geom top)
        },
        "description": "Confined cell, head within bounds → geometric top",
    },
    "cvtdr": {  # Case 3: Convertible, head < bottom (dry cell, tests clamping)
        "top": 10.0,
        "botm": [5.0],
        "strt": [[3.0]],  # head below bottom = dry cell
        "icelltype": [1],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            # Head clamped to bottom: top = max(3.0, 5.0) = 5.0
            # All particles end up at cell bottom (zero saturated thickness)
            0: 5.0,  # bot + 0.0 * (max(head,bot) - bot) = 5.0 + 0.0 * 0 = 5.0
            1: 5.0,  # bot + 0.5 * (max(head,bot) - bot) = 5.0 + 0.5 * 0 = 5.0
            2: 5.0,  # bot + 1.0 * (max(head,bot) - bot) = 5.0 + 1.0 * 0 = 5.0
        },
        "description": "Convertible cell, head < bottom → clamped to bottom",
    },
    "cnflo": {  # Case 4: Confined, head < bottom
        "top": 10.0,
        "botm": [5.0],
        "strt": [[3.0]],  # head below bottom
        "icelltype": [0],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            # Confined cell always uses geometric top, regardless of head
            0: 5.0,  # bot + 0.0 * (top - bot) = 5.0 + 0.0 * 5.0 = 5.0
            1: 7.5,  # bot + 0.5 * (top - bot) = 5.0 + 0.5 * 5.0 = 7.5
            2: 10.0,  # bot + 1.0 * (top - bot) = 5.0 + 1.0 * 5.0 = 10.0
        },
        "description": "Confined cell, head < bottom → geometric top",
    },
    "cvthi": {  # Case 5: Convertible, head > top (tests clamping to top)
        "top": 10.0,
        "botm": [5.0],
        "strt": [[12.0]],  # head above top
        "icelltype": [1],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            # effective top clamped to top: top = min(10.0, 12.0) = 10.0
            0: 5.0,  # bot + 0.0 * (min(top,head) - bot) = 5.0 + 0.0 * 5.0 = 5.0
            1: 7.5,  # bot + 0.5 * (min(top,head) - bot) = 5.0 + 0.5 * 5.0 = 7.5
            2: 10.0,  # bot + 1.0 * (min(top,head) - bot) = 5.0 + 1.0 * 5.0 = 10.0
        },
        "description": "Convertible cell, head > top → clamped to top",
    },
    "cnfhi": {  # Case 6: Confined, head > top
        "top": 10.0,
        "botm": [5.0],
        "strt": [[12.0]],  # head > top
        "icelltype": [0],
        "releasepts": [
            [0, 0, 0, 0, 2.0, 5.0, 0.0],
            [1, 0, 0, 0, 5.0, 5.0, 0.5],
            [2, 0, 0, 0, 8.0, 5.0, 1.0],
        ],
        "expected_z": {
            # Confined cell always uses geometric top, regardless of head
            0: 5.0,  # bot + 0.0 * (top - bot) = 5.0 + 0.0 * 5.0 = 5.0
            1: 7.5,  # bot + 0.5 * (top - bot) = 5.0 + 0.5 * 5.0 = 7.5
            2: 10.0,  # bot + 1.0 * (top - bot) = 5.0 + 1.0 * 5.0 = 10.0
        },
        "description": "Confined cell, head > top → geometric top",
    },
}


def get_case_key(name):
    """Extract the case key from the test name."""
    # Sort keys by length (longest first) to match correctly
    # (e.g., "cvtdr" before "cvt", "cnflo" before "cnf")
    for key in sorted(case_params.keys(), key=len, reverse=True):
        if key in name:
            return key
    return None


def build_gwf_sim(name, ws, mf6):
    """Build the GWF simulation for a given test case."""
    ws = Path(ws)
    gwf_name = get_model_name(name, "gwf")
    case_key = get_case_key(name)
    params = case_params[case_key]

    sim = flopy.mf6.MFSimulation(
        sim_name=gwf_name,
        exe_name=mf6,
        version="mf6",
        sim_ws=ws,
    )

    tdis_rc = [(perlen, nstp, tsmult)]
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    # Use Newton for cases with head below bottom
    use_newton = case_key in ["cvtdr", "cnflo"]
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1e-6,
        inner_dvclose=1e-6,
        linear_acceleration="BICGSTAB" if use_newton else "CG",
    )

    # Use Newton for cases with head below bottom
    use_newton = case_key in ["cvtdr", "cnflo"]
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        save_flows=True,
        newtonoptions="NEWTON" if use_newton else None,
    )

    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=params["top"],
        botm=params["botm"],
    )

    flopy.mf6.ModflowGwfic(gwf, strt=params["strt"])

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_saturation=True,
        save_specific_discharge=True,
        icelltype=params["icelltype"],
        k=1.0,
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwf_name}.cbc",
        head_filerecord=f"{gwf_name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def build_prt_sim(name, gwf_ws, prt_ws, mf6):
    """Build the PRT simulation for a given test case (FMI coupling)."""
    prt_ws = Path(prt_ws)
    gwf_ws = Path(gwf_ws)
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    case_key = get_case_key(name)
    params = case_params[case_key]

    sim = flopy.mf6.MFSimulation(
        sim_name=prt_name,
        exe_name=mf6,
        version="mf6",
        sim_ws=prt_ws,
    )

    tdis_rc = [(perlen, nstp, tsmult)]
    flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

    flopy.mf6.ModflowGwfdis(
        prt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=params["top"],
        botm=params["botm"],
    )

    flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

    prp_track_file = f"{prt_name}.prp.trk"
    prp_track_csv_file = f"{prt_name}.prp.trk.csv"
    flopy.mf6.ModflowPrtprp(
        prt,
        pname="prp1",
        filename=f"{prt_name}_1.prp",
        nreleasepts=len(params["releasepts"]),
        packagedata=params["releasepts"],
        perioddata={0: ["FIRST"]},
        track_filerecord=[prp_track_file],
        trackcsv_filerecord=[prp_track_csv_file],
        local_z=True,
        stop_at_weak_sink=False,
        extend_tracking=True,
    )

    prt_track_file = f"{prt_name}.trk"
    prt_track_csv_file = f"{prt_name}.trk.csv"
    flopy.mf6.ModflowPrtoc(
        prt,
        pname="oc",
        track_filerecord=[prt_track_file],
        trackcsv_filerecord=[prt_track_csv_file],
    )

    gwf_budget_file = gwf_ws / f"{gwf_name}.cbc"
    gwf_head_file = gwf_ws / f"{gwf_name}.hds"
    gwf_grid_file = gwf_ws / f"{gwf_name}.dis.grb"
    flopy.mf6.ModflowPrtfmi(
        prt,
        packagedata=[
            ("GWFGRID", gwf_grid_file),
            ("GWFHEAD", gwf_head_file),
            ("GWFBUDGET", gwf_budget_file),
        ],
    )

    ems = flopy.mf6.ModflowEms(
        sim,
        pname="ems",
        filename=f"{prt_name}.ems",
    )
    sim.register_solution_package(ems, [prt.name])

    return sim


def build_models(idx, test):
    """Build GWF and PRT models for a given test case."""
    name = test.name
    use_exchange = "exg" in name

    if use_exchange:
        # Build GWF and PRT in the same simulation with exchange
        sim = build_gwf_sim(name, test.workspace, test.targets["mf6"])
        gwf_name = get_model_name(name, "gwf")
        prt_name = get_model_name(name, "prt")
        case_key = get_case_key(name)
        params = case_params[case_key]

        # Add PRT model to the same simulation
        prt = flopy.mf6.ModflowPrt(sim, modelname=prt_name)

        flopy.mf6.ModflowGwfdis(
            prt,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=params["top"],
            botm=params["botm"],
        )

        flopy.mf6.ModflowPrtmip(prt, pname="mip", porosity=porosity)

        prp_track_file = f"{prt_name}.prp.trk"
        prp_track_csv_file = f"{prt_name}.prp.trk.csv"
        flopy.mf6.ModflowPrtprp(
            prt,
            pname="prp1",
            filename=f"{prt_name}_1.prp",
            nreleasepts=len(params["releasepts"]),
            packagedata=params["releasepts"],
            perioddata={0: ["FIRST"]},
            track_filerecord=[prp_track_file],
            trackcsv_filerecord=[prp_track_csv_file],
            local_z=True,
            stop_at_weak_sink=False,
            extend_tracking=True,
        )

        prt_track_file = f"{prt_name}.trk"
        prt_track_csv_file = f"{prt_name}.trk.csv"
        flopy.mf6.ModflowPrtoc(
            prt,
            pname="oc",
            track_filerecord=[prt_track_file],
            trackcsv_filerecord=[prt_track_csv_file],
        )

        # Create GWF-PRT exchange
        flopy.mf6.ModflowGwfprt(
            sim,
            exgtype="GWF6-PRT6",
            exgmnamea=gwf_name,
            exgmnameb=prt_name,
            filename=f"{name}.gwfprt",
        )

        # Create EMS package for PRT
        ems = flopy.mf6.ModflowEms(
            sim,
            pname="ems",
            filename=f"{prt_name}.ems",
        )
        sim.register_solution_package(ems, [prt.name])

        return sim, None
    else:
        # Build separate GWF and PRT simulations with FMI
        gwf_sim = build_gwf_sim(name, test.workspace, test.targets["mf6"])
        prt_sim = build_prt_sim(
            name,
            test.workspace,
            test.workspace / "prt",
            test.targets["mf6"],
        )
        return gwf_sim, prt_sim


def check_output(idx, test):
    """Check that particle release z coordinates match expected values."""
    name = test.name
    use_exchange = "exg" in name
    prt_name = get_model_name(name, "prt")
    case_key = get_case_key(name)
    params = case_params[case_key]

    # Determine workspace based on coupling method
    prt_ws = test.workspace if use_exchange else test.workspace / "prt"

    pathlines = pd.read_csv(prt_ws / f"{prt_name}.prp.trk.csv")
    releasepts = pathlines[pathlines["ireason"] == 0]

    print(f"\n{params['description']}")
    for irpt, expected_z_val in params["expected_z"].items():
        actual_z = releasepts[releasepts["irpt"] == irpt + 1]["z"].values[0]
        print(
            f"  Particle {irpt}: "
            f"expected z={expected_z_val:.6f}, "
            f"actual z={actual_z:.6f}"
        )
        assert np.isclose(actual_z, expected_z_val, rtol=1e-5, atol=1e-5), (
            f"Particle {irpt}: expected z={expected_z_val}, got z={actual_z}"
        )


def plot_output(idx, test):
    """Plot particle release points in cross section."""
    name = test.name
    use_exchange = "exg" in name
    gwf_name = get_model_name(name, "gwf")
    prt_name = get_model_name(name, "prt")
    case_key = get_case_key(name)
    params = case_params[case_key]

    # Determine workspace and simulation based on coupling method
    if use_exchange:
        prt_ws = test.workspace
        sim = test.sims[0]
        gwf = sim.get_model(gwf_name)
    else:
        prt_ws = test.workspace / "prt"
        gwf_sim = test.sims[0]
        gwf = gwf_sim.get_model(gwf_name)

    mg = gwf.modelgrid

    prt_track_csv_file = f"{prt_name}.prp.trk.csv"
    track_data = pd.read_csv(prt_ws / prt_track_csv_file, na_filter=False)
    release_points = track_data[track_data["ireason"] == 0]

    # Set up plot
    coupling_method = "exchange" if use_exchange else "FMI"
    use_newton = case_key in ["cvtdr", "cnflo"]
    newton_str = ", Newton" if use_newton else ""
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.set_title(f"{params['description']} ({coupling_method}{newton_str})")

    # Plot grid in cross section
    pmv = flopy.plot.PlotCrossSection(modelgrid=mg, ax=ax, line={"row": 0})
    pmv.plot_grid()

    # Add head line
    head = params["strt"][0][0]
    bot = params["botm"][0]
    top = params["top"]
    icelltype = params["icelltype"][0]

    cell_type_str = "convertible" if icelltype == 1 else "confined"
    if head < bot:
        head_status = "below bottom"
    elif head > top:
        head_status = "above top"
    else:
        head_status = "within bounds"

    ax.axhline(
        y=head,
        color="blue",
        linestyle=":",
        linewidth=2,
        alpha=0.8,
        label=f"Head={head:.1f} ({head_status})",
    )

    # Plot release points
    ax.scatter(
        release_points["x"],
        release_points["z"],
        s=100,
        c="red",
        marker="o",
        label=f"{cell_type_str.capitalize()} cell release points",
        zorder=5,
        edgecolors="black",
    )

    # Annotate points with their local z values
    for _, pt in release_points.iterrows():
        localz = params["releasepts"][int(pt["irpt"]) - 1][6]
        ax.annotate(
            f"lz={localz:.1f}",
            xy=(pt["x"], pt["z"]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7),
        )

    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(prt_ws / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.show()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        targets=targets,
        compare=None,
    )
    test.run()

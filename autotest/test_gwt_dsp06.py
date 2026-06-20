from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["dsp06"]

# Grid and transport controls
nlay = 1
nrow = 61
ncol = 61
delr = 1.0
delc = 1.0
top = 1.0
botm = 0.0

# Flow and transport properties
k = 1.0
porosity = 0.10
hydraulic_gradient = 0.02

# Dispersion parameters
diffc = 0.0
alh = 10.0
alv = 10.0
ath1 = 1.0
ath2 = 1.0
atv = 1.0

# Time discretization
perlen = 80.0
nstp = 80
tsmult = 1.0

# Source concentration
source_concentration = 1.0


def cell_center_xy(i, j):
    """Return x, y coordinates for a structured-grid cell center."""
    x = (j + 0.5) * delr
    y = (nrow - i - 0.5) * delc
    return x, y


def row_col_from_xy(x, y):
    """Return nearest row, col for x, y coordinates."""
    j = round(x / delr - 0.5)
    i = round(nrow - y / delc - 0.5)
    i = int(np.clip(i, 2, nrow - 3))
    j = int(np.clip(j, 2, ncol - 3))
    return i, j


def make_boundary_heads(angle_deg):
    """Create CHD stress-period data for a linear head plane on all boundary cells."""
    theta = np.deg2rad(angle_deg)
    ux = np.cos(theta)
    uy = np.sin(theta)
    h0 = 10.0
    spd = []
    for i in range(nrow):
        for j in range(ncol):
            if i in (0, nrow - 1) or j in (0, ncol - 1):
                x, y = cell_center_xy(i, j)
                h = h0 - hydraulic_gradient * (x * ux + y * uy)
                spd.append(((0, i, j), h))
    return spd


def make_source_cell(angle_deg, upstream_distance=20.0):
    """Place the source upstream of the model center for the selected flow angle."""
    theta = np.deg2rad(angle_deg)
    ux = np.cos(theta)
    uy = np.sin(theta)
    x_center = ncol * delr / 2.0
    y_center = nrow * delc / 2.0
    x_src = x_center - upstream_distance * ux
    y_src = y_center - upstream_distance * uy
    i_src, j_src = row_col_from_xy(x_src, y_src)
    return (0, i_src, j_src)


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name="dsp06", version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=1,
        perioddata=[(perlen, nstp, tsmult)],
    )

    # Create flow models
    gwf_a = flopy.mf6.ModflowGwf(sim, modelname="gwf_dsp06a", save_flows=True)
    gwf_a.set_model_relative_path("gwf_dsp06a")

    gwf_b = flopy.mf6.ModflowGwf(sim, modelname="gwf_dsp06b", save_flows=True)
    gwf_b.set_model_relative_path("gwf_dsp06b")

    # Discretization for both flow models
    for gwf, name in [(gwf_a, "gwf_dsp06a"), (gwf_b, "gwf_dsp06b")]:
        dis = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
            idomain=np.ones((nlay, nrow, ncol), dtype=int),
            filename=f"{name}.dis",
        )
        angle = 60.0 if name == "gwf_dsp06a" else -60.0
        chd_spd = make_boundary_heads(angle)
        strt_arr = np.ones((nlay, nrow, ncol), dtype=float) * 10.0
        for cellid, h in chd_spd:
            strt_arr[cellid] = h

        ic = flopy.mf6.ModflowGwfic(gwf, strt=strt_arr, filename=f"{name}.ic")
        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            icelltype=0,
            k=k,
            save_specific_discharge=True,
            filename=f"{name}.npf",
        )
        chd = flopy.mf6.ModflowGwfchd(
            gwf,
            stress_period_data={0: chd_spd},
            save_flows=False,
            pname="CHD-1",
            filename=f"{name}.chd",
        )

    # Output control for both flow models
    for gwf, name in [(gwf_a, "gwf_dsp06a"), (gwf_b, "gwf_dsp06b")]:
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{name}.cbc",
            head_filerecord=f"{name}.hds",
            saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
            filename=f"{name}.oc",
        )

    # Create iterative model solution for flow and register both flow models
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1e-10,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=100,
        inner_dvclose=1e-10,
        rcloserecord=1e-10,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=1.0,
        filename="flow.ims",
    )
    sim.register_ims_package(imsgwf, [gwf_a.name, gwf_b.name])

    # Create transport models
    gwt_a = flopy.mf6.MFModel(sim, model_type="gwt6", modelname="gwt_dsp06a")
    gwt_a.name_file.save_flows = True
    gwt_a.set_model_relative_path("gwt_dsp06a")

    gwt_b = flopy.mf6.MFModel(sim, model_type="gwt6", modelname="gwt_dsp06b")
    gwt_b.name_file.save_flows = True
    gwt_b.set_model_relative_path("gwt_dsp06b")

    # Discretization, IC, ADV, DSP, MST, SSM, OC for both transport models
    for gwt, name in [(gwt_a, "gwt_dsp06a"), (gwt_b, "gwt_dsp06b")]:
        dis = flopy.mf6.ModflowGwtdis(
            gwt,
            nlay=nlay,
            nrow=nrow,
            ncol=ncol,
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
            idomain=1,
            filename=f"{name}.dis",
        )
        ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{name}.ic")
        adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{name}.adv")
        dsp = flopy.mf6.ModflowGwtdsp(
            gwt,
            xt3d_off=False,
            diffc=diffc,
            alh=alh,
            alv=alv,
            ath1=ath1,
            ath2=ath2,
            atv=atv,
            filename=f"{name}.dsp",
        )
        mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{name}.mst")
        ssm = flopy.mf6.ModflowGwtssm(gwt, sources=[[]], filename=f"{name}.ssm")
        oc = flopy.mf6.ModflowGwtoc(
            gwt,
            budget_filerecord=f"{name}.cbc",
            concentration_filerecord=f"{name}.ucn",
            saverecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
            printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
            filename=f"{name}.oc",
        )

    # Constant concentration sources
    source_cell_a = make_source_cell(60.0)
    cnc_a = flopy.mf6.ModflowGwtcnc(
        gwt_a,
        stress_period_data={0: [[source_cell_a, source_concentration]]},
        save_flows=False,
        pname="CNC-1",
        filename="gwt_dsp06a.cnc",
    )
    source_cell_b = make_source_cell(-60.0)
    cnc_b = flopy.mf6.ModflowGwtcnc(
        gwt_b,
        stress_period_data={0: [[source_cell_b, source_concentration]]},
        save_flows=False,
        pname="CNC-1",
        filename="gwt_dsp06b.cnc",
    )

    # Create iterative model solution for transport and register both transport models
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        outer_dvclose=1e-10,
        outer_maximum=100,
        under_relaxation="NONE",
        inner_maximum=300,
        inner_dvclose=1e-10,
        rcloserecord=1e-10,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=1.0,
        filename="transport.ims",
    )
    sim.register_ims_package(imsgwt, [gwt_a.name, gwt_b.name])

    # GWF-GWT exchanges
    gwfgwt_a = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="gwf_dsp06a",
        exgmnameb="gwt_dsp06a",
        filename="gwf_dsp06a_gwt_dsp06a.gwfgwt",
    )
    gwfgwt_b = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea="gwf_dsp06b",
        exgmnameb="gwt_dsp06b",
        filename="gwf_dsp06b_gwt_dsp06b.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    # Define paths to the concentration files
    ws = Path(test.workspace)

    fpth_a = ws / "gwt_dsp06a" / "gwt_dsp06a.ucn"
    if not fpth_a.exists():
        fpth_a = ws / "gwt_dsp06a.ucn"

    fpth_b = ws / "gwt_dsp06b" / "gwt_dsp06b.ucn"
    if not fpth_b.exists():
        fpth_b = ws / "gwt_dsp06b.ucn"

    # Load concentrations
    cobj_a = flopy.utils.HeadFile(fpth_a, precision="double", text="CONCENTRATION")
    c_pos = cobj_a.get_alldata()[-1]  # Final time step, shape (nlay, nrow, ncol)

    cobj_b = flopy.utils.HeadFile(fpth_b, precision="double", text="CONCENTRATION")
    c_neg = cobj_b.get_alldata()[-1]  # Final time step, shape (nlay, nrow, ncol)

    # Diagnostic checks to ensure the plume has developed and spread
    assert c_pos.max() > 0.01, f"Case A plume too weak: max={c_pos.max()}"
    assert c_neg.max() > 0.01, f"Case B plume too weak: max={c_neg.max()}"
    assert np.count_nonzero(c_pos > 0.001) > 10, "Case A plume too small"
    assert np.count_nonzero(c_neg > 0.001) > 10, "Case B plume too small"

    # Flip the negative flow concentration along the row axis (axis 1 of 3D array)
    c_neg_flipped = np.flip(c_neg, axis=1)

    # Compare
    diff = np.abs(c_pos - c_neg_flipped)
    max_diff = np.max(diff)

    print(f"max_diff with fixed code: {max_diff}")

    assert max_diff < 1e-7, (
        f"DSP horizontal angle bug detected! Max plume difference: {max_diff}"
    )


def plot_output(idx, test):
    print("Plotting results")
    import matplotlib.pyplot as plt

    ws = Path(test.workspace)

    fpth_a = ws / "gwt_dsp06a" / "gwt_dsp06a.ucn"
    if not fpth_a.exists():
        fpth_a = ws / "gwt_dsp06a.ucn"

    fpth_b = ws / "gwt_dsp06b" / "gwt_dsp06b.ucn"
    if not fpth_b.exists():
        fpth_b = ws / "gwt_dsp06b.ucn"

    # Load concentrations
    cobj_a = flopy.utils.HeadFile(fpth_a, precision="double", text="CONCENTRATION")
    c_pos = cobj_a.get_alldata()[-1][
        0
    ]  # Final time step, first layer, shape (nrow, ncol)

    cobj_b = flopy.utils.HeadFile(fpth_b, precision="double", text="CONCENTRATION")
    c_neg = cobj_b.get_alldata()[-1][
        0
    ]  # Final time step, first layer, shape (nrow, ncol)

    c_neg_flipped = np.flip(c_neg, axis=0)  # flip row axis (up-down)
    diff = np.abs(c_pos - c_neg_flipped)

    print(c_pos.shape, c_neg.shape, c_neg_flipped.shape, diff.shape)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    # Plot GWT1 (60 degrees)
    im0 = axes[0].imshow(
        np.flipud(c_pos),
        origin="lower",
        extent=(0, ncol * delr, 0, nrow * delc),
        vmin=0.0,
        vmax=1.0,
    )
    # Source cell for GWT1
    _, i_src_a, j_src_a = make_source_cell(60.0)
    x_src_a, y_src_a = cell_center_xy(i_src_a, j_src_a)
    axes[0].plot(x_src_a, y_src_a, marker="o", markersize=6, color="red")
    # Arrow for GWT1
    theta_a = np.deg2rad(60.0)
    axes[0].arrow(
        ncol * delr * 0.80,
        nrow * delc * 0.15,
        8.0 * np.cos(theta_a),
        8.0 * np.sin(theta_a),
        width=0.3,
        length_includes_head=True,
        color="black",
    )
    axes[0].set_title("GWT1 (60°)")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    # Plot GWT2 (-60 degrees)
    im1 = axes[1].imshow(
        np.flipud(c_neg),
        origin="lower",
        extent=(0, ncol * delr, 0, nrow * delc),
        vmin=0.0,
        vmax=1.0,
    )
    # Source cell for GWT2
    _, i_src_b, j_src_b = make_source_cell(-60.0)
    x_src_b, y_src_b = cell_center_xy(i_src_b, j_src_b)
    axes[1].plot(x_src_b, y_src_b, marker="o", markersize=6, color="red")
    # Arrow for GWT2
    theta_b = np.deg2rad(-60.0)
    axes[1].arrow(
        ncol * delr * 0.80,
        nrow * delc * 0.15,
        8.0 * np.cos(theta_b),
        8.0 * np.sin(theta_b),
        width=0.3,
        length_includes_head=True,
        color="black",
    )
    axes[1].set_title("GWT2 (-60°)")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    # Plot Difference (GWT1 - flipped GWT2)
    im2 = axes[2].imshow(
        np.flipud(diff),
        origin="lower",
        extent=(0, ncol * delr, 0, nrow * delc),
    )
    axes[2].set_title("Absolute Difference (GWT1 - flipped GWT2)")
    axes[2].set_xlabel("x")
    axes[2].set_ylabel("y")
    axes[2].set_aspect("equal")

    fig.colorbar(im1, ax=axes[:2], label="Concentration")
    fig.colorbar(im2, ax=axes[2], label="Difference")

    fname = test.workspace / "gwt_dsp06.png"
    print(f"Saving plot to {fname}")
    plt.savefig(fname)
    plt.close()


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
    )
    test.run()

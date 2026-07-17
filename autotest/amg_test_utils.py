"""Shared helpers for the IMS AMG preconditioner tests.

Collects the pieces that were duplicated across the ``test_*_ims_amg*.py``
files: pandas-based readers for the IMS iteration CSV output, and the shared
heterogeneous confined GWF benchmark (geometry, checkerboard K field, and IMS
keyword arguments).
"""

import flopy
import numpy as np
import pandas as pd

# geometry of the shared heterogeneous confined GWF benchmark
NLAY, NROW, NCOL = 3, 10, 20
DELR = DELC = 100.0
TOP = 0.0
BOTM = [-5.0, -10.0, -20.0]
STRT = 5.0
CHD_LEFT = 10.0
CHD_RIGHT = 5.0


def count_inner_iters(csv_path):
    """Total inner (linear) iterations from an IMS CSV output file.

    Works with either the inner or outer IMS CSV; ``total_inner_iterations`` is
    a cumulative counter, so its maximum is the run total.
    """
    return int(pd.read_csv(csv_path)["total_inner_iterations"].max())


def inner_iters_per_timestep(csv_path):
    """Inner iterations summed per (kper, kstp) from an IMS outer CSV.

    Lets callers compare convergence time step by time step rather than as a
    single simulation total.
    """
    df = pd.read_csv(csv_path)
    per = df.groupby(["kper", "kstp"])["inner_iterations"].sum().sort_index()
    return per.to_numpy()


def build_checker_gwf(sim, name):
    """Add the shared 3-layer confined GWF model to *sim*.

    Layer base K values span three orders of magnitude and a checkerboard
    pattern within each layer puts small-K cells directly next to large-K cells
    (100:1 lateral contrast). CHD boundaries are set on the left and right faces
    of every layer.
    """
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=NLAY,
        nrow=NROW,
        ncol=NCOL,
        delr=DELR,
        delc=DELC,
        top=TOP,
        botm=BOTM,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=STRT)
    k_layer = np.array([10.0, 0.1, 1.0])  # confining unit 100x lower than layer 1
    rows = np.arange(NROW)
    cols = np.arange(NCOL)
    checker = np.where((rows[:, None] + cols[None, :]) % 2 == 0, 10.0, 0.1)
    k = k_layer[:, None, None] * checker[None, :, :]
    flopy.mf6.ModflowGwfnpf(gwf, k=k, k33=k * 0.1)
    spd = [[(il, ir, 0), CHD_LEFT] for il in range(NLAY) for ir in range(NROW)]
    spd += [[(il, ir, NCOL - 1), CHD_RIGHT] for il in range(NLAY) for ir in range(NROW)]
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=spd)
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )
    return gwf


def base_ims_kwargs(accel, relax, inner_maximum=500, preconditioner_levels=5):
    """IMS keyword arguments shared by the AMG/ILU comparison tests.

    relax = 1.0 gives a fixed (undamped) AMG cycle or MILU/MILUT for ILU;
    relax = 0.0 gives adaptive AMG damping or ILU0/ILUT for ILU.
    """
    return dict(
        print_option="ALL",
        outer_dvclose=1.0e-9,
        outer_maximum=100,
        inner_dvclose=1.0e-12,
        inner_maximum=inner_maximum,
        rcloserecord="1.0e-9",
        linear_acceleration=accel,
        preconditioner_levels=preconditioner_levels,
        relaxation_factor=relax,
    )

"""
Compare a single anisotropic model with the same domain split in two
and joined by a GWF-GWF exchange. Cells are convertible with the water
table inside the layer, so the exchange conductance is recomputed every
outer iteration from the effective K in the direction given by ANGLDEGX.

Case "x" splits the domain in x (ANGLDEGX = 0), case "y" splits it in y
(ANGLDEGX = 90). Heads in the two-model simulation match the reference.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["gwfgwf_angx", "gwfgwf_angy"]
splits = ["x", "y"]

nrow, ncol = 10, 10
delr = delc = 100.0
top, bot = 10.0, 0.0
k11, k22 = 1.0, 0.1
h_nw, h_se = 8.0, 2.0
hclose = 1e-9


def add_model(sim, name, nrow, ncol, chd, yorigin=0.0, xorigin=0.0):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=bot,
        xorigin=xorigin,
        yorigin=yorigin,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, k=k11, k22=k22)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chd)
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "ALL")]
    )


def build_models(idx, test):
    split = splits[idx]
    sim = flopy.mf6.MFSimulation(sim_name=cases[idx], sim_ws=test.workspace)
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=hclose,
        outer_maximum=200,
        inner_maximum=300,
        inner_dvclose=hclose,
        rcloserecord=1e-6,
        linear_acceleration="BICGSTAB",
    )

    add_model(
        sim, "ref", nrow, ncol, [[(0, 0, 0), h_nw], [(0, nrow - 1, ncol - 1), h_se]]
    )

    if split == "x":
        # model a is west, model b is east; connections face +x
        nc = ncol // 2
        add_model(sim, "a", nrow, nc, [[(0, 0, 0), h_nw]])
        add_model(
            sim, "b", nrow, nc, [[(0, nrow - 1, nc - 1), h_se]], xorigin=nc * delr
        )
        exgdata = [
            [(0, i, nc - 1), (0, i, 0), 1, delr / 2, delr / 2, delc, 0.0]
            for i in range(nrow)
        ]
    else:
        # model a is south, model b is north; connections face +y
        nr = nrow // 2
        add_model(sim, "a", nr, ncol, [[(0, nr - 1, ncol - 1), h_se]])
        add_model(sim, "b", nr, ncol, [[(0, 0, 0), h_nw]], yorigin=nr * delc)
        exgdata = [
            [(0, 0, j), (0, nr - 1, j), 1, delc / 2, delc / 2, delr, 90.0]
            for j in range(ncol)
        ]

    flopy.mf6.ModflowGwfgwf(
        sim,
        exgtype="GWF6-GWF6",
        exgmnamea="a",
        exgmnameb="b",
        nexg=len(exgdata),
        exchangedata=exgdata,
        auxiliary=["ANGLDEGX"],
    )
    return sim, None


def check_output(idx, test):
    ws = test.workspace
    ref = flopy.utils.HeadFile(ws / "ref.hds").get_data()
    a = flopy.utils.HeadFile(ws / "a.hds").get_data()
    b = flopy.utils.HeadFile(ws / "b.hds").get_data()
    if splits[idx] == "x":
        split = np.concatenate([a, b], axis=2)
    else:
        split = np.concatenate([b, a], axis=1)
    assert (ref > bot).all() and (ref < top).all(), "water table should be in layer"
    assert np.allclose(ref, split, atol=1e-6), (
        f"max head diff {np.abs(ref - split).max()}"
    )


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

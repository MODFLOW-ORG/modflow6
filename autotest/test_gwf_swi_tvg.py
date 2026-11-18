"""
Test the Time-Varying Grid (TVG) Utility as used by the SWI Package.

The saltwater head can be supplied two ways: as the scalar SALTWATER_HEAD
option, or as a full-grid array through a TVG6 file. For a grid whose cells
all carry the same value the two must give identical results, which checks
that TVG reaches every node rather than only the first layer. The grid is
three layers deliberately: a first-layer-only fill fails this test.

A second case gives a different value in each layer, through the LAYERED
keyword, and confirms the interface responds layer by layer.
"""

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["swi-tvg-uniform", "swi-tvg-layered"]

nlay, nrow, ncol = 3, 1, 6
delr = delc = 100.0
top = 0.0
botm = [-10.0, -20.0, -30.0]
k = 1.0
recharge = 1.0e-4
hsalt = 0.05
hsalt_layered = [0.05, 0.10, 0.15]


def _base(ws, name):
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, nper=1)
    flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        linear_acceleration="bicgstab",
        outer_maximum=200,
    )
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname="mymodel", save_flows=True, newtonoptions="newton"
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=k)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 0.0]])
    flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
    return sim, gwf


def build_models(idx, test):
    ws = pl.Path(test.workspace)

    # reference: the same value supplied through the scalar option
    if idx == 0:
        ref, gwfr = _base(str(ws / "ref"), "ref")
        flopy.mf6.ModflowGwfswi(
            gwfr, zeta_filerecord="mymodel.zta", saltwater_head=hsalt
        )
        flopy.mf6.ModflowGwfoc(
            gwfr, head_filerecord="mymodel.hds", saverecord=[("HEAD", "ALL")]
        )
        ref.write_simulation(silent=True)

    sim, gwf = _base(str(ws), "mymodel")
    swi = flopy.mf6.ModflowGwfswi(
        gwf, zeta_filerecord="mymodel.zta", tvg6_filerecord="mymodel.swi.tvg"
    )
    if idx == 0:
        arr = np.full((nlay, nrow, ncol), hsalt)
    else:
        arr = np.empty((nlay, nrow, ncol))
        for kk in range(nlay):
            arr[kk] = hsalt_layered[kk]
    swi.tvg.initialize(auxiliary=["saltwater_head"], aux={0: [arr]})
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord="mymodel.hds", saverecord=[("HEAD", "ALL")]
    )
    return sim, None


def _zeta(ws):
    return (
        flopy.utils.HeadFile(pl.Path(ws) / "mymodel.zta", text="zeta")
        .get_data()
        .flatten()
    )


def check_output(idx, test):
    ws = pl.Path(test.workspace)
    gwf = test.sims[0].gwf[0]
    head = gwf.output.head().get_data().flatten()
    zeta = _zeta(ws)

    if idx == 0:
        # run the scalar-option reference and require an identical answer
        ref = flopy.mf6.MFSimulation.load(
            sim_ws=str(ws / "ref"), exe_name=test.targets["mf6"], verbosity_level=0
        )
        success, buff = ref.run_simulation(silent=True)
        assert success, "reference simulation failed"
        href = ref.get_model("mymodel").output.head().get_data().flatten()
        zref = _zeta(ws / "ref")
        assert np.allclose(head, href, atol=1.0e-12), (
            "TVG head differs from the scalar SALTWATER_HEAD option: "
            f"max |diff| = {np.abs(head - href).max():.3e}"
        )
        assert np.allclose(zeta, zref, atol=1.0e-12), (
            "TVG zeta differs from the scalar SALTWATER_HEAD option, which "
            "means TVG did not fill every node: "
            f"max |diff| = {np.abs(zeta - zref).max():.3e}"
        )
        print(f"TVG matches the scalar option to {np.abs(zeta - zref).max():.1e}")
    else:
        # a different value per layer must move the interface per layer
        means = [zeta[kk * ncol : (kk + 1) * ncol].mean() for kk in range(nlay)]
        assert len(set(np.round(means, 6))) == nlay, (
            f"layered saltwater head did not give a distinct interface per "
            f"layer: layer means {means}"
        )
        print(f"layered zeta means: {[round(float(m), 4) for m in means]}")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
    )
    test.run()

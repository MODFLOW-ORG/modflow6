"""
Test that an embedded lake's table sets the bottom of the embedded lake's
own connection, not another lake's.

An embedded lake has a single connection, whose bottom comes from the
lake's table instead of the cell bottom. LAK mistakenly looked up the
connection by lake number instead of connection number. The two are the
same only when every lake before the embedded lake has one connection,
so otherwise the table bottom went to another lake's connection, and the
embedded lake's kept the cell bottom.

Both models have a lake with two identical horizontal connections and an
embedded lake, in opposite orders. In hfirst the embedded lake is lake 2
but connection 3, so it was affected. In efirst it is lake 1 and
connection 1, so it wasn't. The conductances and flows must be the same
in both models, and the horizontal lake's two connections must have the
same conductance. The embedded cell's head is below the table bottom, so
the embedded lake's flow depends on its bottom too.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["lak-embedded-belev"]
models = {"hfirst": ("horizontal", "embedded"), "efirst": ("embedded", "horizontal")}

top = 10.0
bot = -10.0
horizontal_stage = 8.0
embedded_stage = 6.0
# stage, volume, sarea, warea; the empty first row puts the lake
# bottom at 2.0, above the cell bottom (-10.0)
embedded_table = [
    [2.0, 0.0, 0.0, 0.0],
    [4.0, 50.0, 50.0, 50.0],
    [10.0, 400.0, 100.0, 100.0],
]


def add_model(sim, name, order):
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=4, delr=100.0, delc=100.0, top=top, botm=bot
    )
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=10.0)
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data=[[(0, 0, 0), 0.0]])

    packagedata = []
    connectiondata = []
    perioddata = []
    obs = []
    for ilak, kind in enumerate(order):
        if kind == "horizontal":
            stage = horizontal_stage
            conns = [
                [(0, 0, j), "horizontal", 0.1, 0.0, 5.0, 10.0, 10.0] for j in (1, 2)
            ]
        else:
            stage = embedded_stage
            conns = [[(0, 0, 3), "embeddedh", 0.1, 0.0, 0.0, 10.0, 0.0]]
        packagedata.append([ilak, stage, len(conns)])
        for iconn, conn in enumerate(conns):
            connectiondata.append([ilak, iconn, *conn])
            obs += [
                (f"{kind[0]}c{iconn + 1}", "conductance", ilak + 1, iconn + 1),
                (f"{kind[0]}q{iconn + 1}", "lak", ilak + 1, iconn + 1),
            ]
        perioddata.append([ilak, "status", "constant"])
    ilak = order.index("embedded")
    lak = flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=len(order),
        ntables=1,
        tables=[[ilak, f"{name}.lak.tab"]],
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata={0: perioddata},
        observations={f"{name}.lak.obs.csv": obs},
    )
    flopy.mf6.ModflowUtllaktab(
        gwf,
        nrow=len(embedded_table),
        ncol=len(embedded_table[0]),
        table=embedded_table,
        filename=f"{name}.lak.tab",
        pname="tab",
        parent_file=lak,
    )
    return gwf


def build_models(idx, test):
    name = cases[idx]
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim)
    flopy.mf6.ModflowIms(sim)
    for model, order in models.items():
        add_model(sim, model, order)
    return sim, None


def check_output(idx, test):
    obs = {
        model: np.genfromtxt(
            test.workspace / f"{model}.lak.obs.csv", delimiter=",", names=True
        )
        for model in models
    }
    for col in ("HC1", "HC2", "HQ1", "HQ2", "EC1", "EQ1"):
        a, b = obs["hfirst"][col], obs["efirst"][col]
        assert np.allclose(a, b), f"{col}: hfirst {a} != efirst {b}"
    for model, values in obs.items():
        c1, c2 = values["HC1"], values["HC2"]
        assert np.allclose(c1, c2), (
            f"{model}: horizontal lake conductances differ: {c1} != {c2}"
        )


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

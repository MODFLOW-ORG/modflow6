# Test evap in SFR reaches (no interaction with gwf)


import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "sfr-bedk",
    "sfr-bedkts",
    "sfr-rhkts",
]


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    length_units = "m"
    time_units = "sec"

    nrow = 1
    ncol = 1
    nlay = 1
    delr = delc = 1.0

    nper = 4
    tds_spd = [
        (1.0, 1, 1.0),
        (1.0, 1, 1.0),
        (1.0, 1, 1.0),
        (1.0, 1, 1.0),
    ]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units=time_units,
        nper=nper,
        perioddata=tds_spd,
    )
    flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
    )
    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=-100.0,
    )
    flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,  # >0 means saturated thickness varies with computed head
    )
    flopy.mf6.ModflowGwfic(gwf, strt=1.0)
    flopy.mf6.ModflowGwfghb(gwf, stress_period_data=[((0, 0, 0), 1.0, 1e6)])

    # sfr data
    nreaches = 2

    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #        <man> <ncon> <ustrf> <ndv>
    if idx < 2:
        rhk1 = 0.0
        rhk2 = 0.0
    else:
        rhk1 = "bedk1"
        rhk2 = "bedk2"
    package_data = [
        (0, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, rhk1, 0.001, 0, 0.0, 0),
        (1, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, rhk2, 0.001, 0, 0.0, 0),
    ]
    connection_data = [
        (0,),
        (1,),
    ]

    if idx == 0:
        timeseries = False
        sfr_spd = {
            0: [
                (1, "bedk", 10.0),
            ],
            1: [
                (0, "bedk", 1.0),
                (1, "bedk", 5.0),
            ],
            2: [
                (0, "bedk", 5.0),
                (1, "bedk", 1.0),
            ],
            3: [
                (0, "bedk", 10.0),
                (1, "bedk", 0.0),
            ],
        }
    else:
        timeseries = True
        ts_names = ["bedk1", "bedk2"]
        ts_methods = ["linearend"] * len(ts_names)
        ts_data = [
            (0.0, 0.0, 10.0),
            (1.0, 0.0, 10.0),
            (2.0, 1.0, 5.0),
            (3.0, 5.0, 1.0),
            (4.0, 10.0, 0.0),
        ]
        if idx < 2:
            sfr_spd = {
                0: [
                    (0, "bedk", "bedk1"),
                    (1, "bedk", "bedk2"),
                ]
            }
        else:
            sfr_spd = None

    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwfr1", "sfr", (0,)),
            ("gwfr2", "sfr", (1,)),
        ],
        "filename": f"{name}.sfr.obs",
    }

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=1.0,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=sfr_spd,
        observations=sfr_obs,
        pname="SFR-1",
    )

    if timeseries:
        sfr.ts.initialize(
            filename=f"{name}.sfr.ts",
            timeseries=ts_data,
            time_series_namerecord=ts_names,
            interpolation_methodrecord=ts_methods,
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("head", "all"), ("budget", "all")],
    )

    return sim, None


def check_output(idx, test):
    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwf = sim.get_model()
    sfr = gwf.get_package("SFR-1")
    obs_data = sfr.output.obs().get_data()
    o1 = obs_data["GWFR1"]
    o2 = obs_data["GWFR2"][::-1]
    assert np.allclose(o1, o2), f"GWFR1 ({o1}) not equal to reversed GWFR2 ({o2})"


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

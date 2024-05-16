# Test groundwater discharge to a stream and then go on to test
# that a transport model with a single reach works.

import math
import pathlib as pl

import flopy
import numpy as np
import pytest

from framework import TestFramework

cases = ["sfr-gwfout", "sfr-gwf-trnsprt"]

length_units = "m"
time_units = "sec"

nrow = 1
ncol = 1
nlay = 1
delr = delc = 1.0

nouter, ninner = 100, 300
dvclose, rclose, relax = 1e-10, 1e-10, 1.0


def add_gwt_model(sim, gwtname):
    gwt = flopy.mf6.ModflowGwt(
        sim, modelname=gwtname, model_nam_file="{}.nam".format(gwtname)
    )
    gwt.name_file.save_flows = True

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=dvclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=dvclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwtname),
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=-100.0,
    )

    flopy.mf6.ModflowGwtic(gwt, strt=1.0)
    flopy.mf6.ModflowGwtmst(gwt, porosity=0.2, filename=f"{gwtname}.mst")

    sourcerecarray = [("GHB", "AUX", "CONCENTRATION")]
    flopy.mf6.ModflowGwtssm(
        gwt,
        sources=sourcerecarray,
        pname="SSM",
        filename="{}.ssm".format(gwtname),
    )

    # Instantiate Streamflow Transport package
    sft_packagedata = []
    t = (0, 1.0)
    sft_packagedata.append(t)

    sft_perioddata = []
    sft_perioddata.append((0, "INFLOW", 0.0))
    flwpckname = "SFR"

    flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sft_packagedata,
        reachperioddata=sft_perioddata,
        flow_package_name=flwpckname,
        pname="SFT",
        filename="{}.sft".format(gwtname),
    )

    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    return sim


def build_models(idx, test):
    # Base simulation and model name and workspace
    ws = test.workspace
    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=ws, exe_name="mf6", version="mf6"
    )
    flopy.mf6.ModflowTdis(sim, time_units=time_units)
    flopy.mf6.ModflowIms(
        sim,
        inner_dvclose=1e-5,
        inner_hclose=1e-6,
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
    flopy.mf6.ModflowGwfghb(
        gwf,
        auxiliary="CONCENTRATION",
        stress_period_data=[((0, 0, 0), 1.0, 1e6, 1.0)],
        pname="GHB",
    )

    # sfr data
    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> <man> <ncon> <ustrf> <ndv>
    package_data = [
        (0, (0, 0, 0), delr, 1.0, 1e-3, 0.0, 1.0, 1.0, 0.001, 0, 0.0, 0)
    ]
    connection_data = [(0)]

    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwf", "sfr", (0,)),
            ("outflow", "ext-outflow", (0,)),
            ("depth", "depth", (0,)),
        ],
        "filename": name + ".sfr.obs",
    }

    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        length_conversion=1.0,
        time_conversion=1.0,
        nreaches=1,
        packagedata=package_data,
        connectiondata=connection_data,
        observations=sfr_obs,
        pname="SFR",
    )

    if idx > 0:
        gwtname = "gwt-sft"
        sim = add_gwt_model(sim, gwtname)

        # Add the flow-transport exchanges
        flopy.mf6.ModflowGwfgwt(
            sim,
            exgtype="GWF6-GWT6",
            exgmnamea=name,
            exgmnameb=gwtname,
            pname="GWFGWT1",
            filename="{}.gwfgwt1".format(gwtname),
        )

    return sim, None


def check_output(idx, test):
    answer = np.array(
        [
            1.0,
            -0.92094535738673577,
            -0.92094535738673577,
            0.79053721667952215e-1,
        ]
    )
    obs_pth = pl.Path(f"{test.workspace}/{cases[idx]}.sfr.obs.csv")
    sim_data = flopy.utils.Mf6Obs(obs_pth).get_data()
    data_names = sim_data.dtype.names
    for ct, name in enumerate(data_names):
        assert np.allclose(
            sim_data[name][0], answer[ct]
        ), f"simulated sfr {name} results do not match answer"

    if idx > 0:
        sft_ans = np.array([[[[1.0]]]])

        gwtname = "gwt-sft"
        sft_obs_fl = pl.Path(f"{test.workspace}/{gwtname}.sft.bin")
        try:
            # load simulated concentration in SFT
            cobj = flopy.utils.HeadFile(
                sft_obs_fl, text="CONCENTRATION"  # precision="double"
            )
            sim_conc_sft = cobj.get_alldata()
        except:
            assert (
                False
            ), f'could not load concentration data from "{sft_obs_fl}"'

        gwt_sim_conc = pl.Path(f"{test.workspace}/{gwtname}.ucn")
        try:
            # load simulated concentration of groundwater
            cobj = flopy.utils.HeadFile(
                gwt_sim_conc, text="CONCENTRATION"  # precision="double"
            )
            conc_gw = cobj.get_alldata()
        except:
            assert (
                False
            ), f'could not load temperature data from "{sft_obs_fl}"'

        msg0 = "The simulation is not matching the established answer"
        msg1 = (
            "Groundwater discharge is the only source of flow in the "
            "channel. Thus, the gw and sft water should have the same "
            "concentration, but don't"
        )
        assert np.allclose(sft_ans, sim_conc_sft), msg0
        assert np.allclose(conc_gw, sim_conc_sft), msg1


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

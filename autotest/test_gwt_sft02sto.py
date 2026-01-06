"""
Simple one-layer model with sfr on top. Purpose is to test transport in a one-d
sfr network with corrected storage term. SFR inflows and rainfall in the sfr network
were specified to result in a defined stage. Initially the inflow to the reaches is
only from the upstream boundary. Later in the simulation rainfall is the primary
source of inflow to the reach. There is no flow between  the stream and the aquifer.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sft_02sto"]

RHK = 0.0
MANNINGS = 0.01
RLEN = 1.0
RWID = 1.0
RAREA = RLEN * RWID
SLOPE = 0.001
INFLOW_CONC = 1.0

SFR_DEPTH = np.array(
    [
        0.1,
        1.0,
        2.0,
        3.0,
        4.0,
        3.0,
        2.0,
        1.0,
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
        3.0,
        2.0,
        1.0,
        0.1,
    ]
)
INFLOW_DEPTH = np.array(
    [
        0.1,
        1.0,
        2.0,
        3.0,
        4.0,
        3.0,
        2.0,
        1.0,
        0.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        0.1,
    ]
)

INFLOW_TOTAL = (1.0 / MANNINGS) * RWID * SFR_DEPTH ** (5 / 3) * np.sqrt(SLOPE)
INFLOW_RATE = (1.0 / MANNINGS) * RWID * INFLOW_DEPTH ** (5 / 3) * np.sqrt(SLOPE)
RECHARGE_RATE = INFLOW_TOTAL - INFLOW_RATE
RECHARGE_DEPTH = RECHARGE_RATE / RAREA
Q_MP = INFLOW_RATE + 0.5 * RECHARGE_RATE
SFR_DEPTH_MP = (Q_MP * MANNINGS / (RWID * np.sqrt(SLOPE))) ** (3 / 5)
dtype = [("D01", float), ("D02", float)]
dtype = [("V01", float), ("V02", float)]
SFR_VOLUME_TS = np.array(
    [(d0 * RAREA, d1 * RAREA) for d0, d1 in zip(SFR_DEPTH_MP, SFR_DEPTH)], dtype=dtype
)


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]

    nlay, nrow, ncol = 1, 1, 1
    delr = RLEN
    delc = 1.0
    top = 0.0
    botm = [-100.0]
    nper = len(INFLOW_TOTAL)
    perioddata = [(1.0, 1, 1.0)] * nper

    # <ifno> <cellid> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> <man> <ncon> <ustrf> <ndv>
    conn_data = [[0, -1], [1, 0]]
    package_data = [
        (0, (-1, -1, -1), RLEN, RWID, SLOPE, top, 1.0, RHK, MANNINGS, 1, 1.0, 0),
        (1, (-1, -1, -1), RLEN, RWID, SLOPE, top, 1.0, RHK, MANNINGS, 1, 1.0, 0),
    ]
    nreaches = len(package_data)

    sfr_dict = {}
    for n in range(nper):
        sfr_spd = [
            (0, "inflow", INFLOW_RATE[n]),
            (0, "rainfall", RECHARGE_DEPTH[n]),
        ]
        sfr_dict[n] = sfr_spd

    # build the simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=ws,
    )
    tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)

    # gwf model
    gwf_name = f"{name}_gwf"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwf_name,
        model_nam_file=f"{gwf_name}.nam",
        newtonoptions="newton under_relaxation",
    )

    ims_gwf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
        linear_acceleration="bicgstab",
        outer_maximum=200,
        inner_maximum=100,
        pname="ims-gwf",
        filename=f"{gwf_name}.ims",
    )

    sim.register_ims_package(ims_gwf, [gwf_name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwf_name}.dis",
    )

    npf = flopy.mf6.ModflowGwfnpf(gwf, k=10.0, icelltype=1)

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-5,
        sy=0.25,
        transient={0: True},
        filename=f"{gwf_name}.sto",
    )

    ic = flopy.mf6.ModflowGwfic(
        gwf,
        strt=top,
        filename=f"{gwf_name}.ic",
    )

    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data={0: [[(0, 0, 0), 0]]},
        filename=f"{gwf_name}.chd",
    )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("budget", "all"), ("head", "all")],
        filename=f"{gwf_name}.oc",
    )

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        time_conversion=1.0,
        print_input=True,
        print_stage=True,
        print_flows=True,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=conn_data,
        perioddata=sfr_dict,
        pname="SFR-1",
        filename=f"{gwf_name}.sfr",
    )

    sfr_obs = [(f"DEPTH{n + 1:02d}", "DEPTH", (n,)) for n in range(nreaches)]
    sfr_obs = {f"{gwf_name}.sfr.obs.csv": sfr_obs}
    sfr.obs.initialize(continuous=sfr_obs)

    # gwt model
    gwt_name = f"{name}_gwt"
    gwt = flopy.mf6.ModflowGwt(
        sim, modelname=gwt_name, model_nam_file=f"{gwt_name}.nam", save_flows=True
    )

    ims_gwt = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
        linear_acceleration="bicgstab",
        outer_maximum=200,
        inner_maximum=100,
        filename=f"{gwt_name}.ims",
        pname="ims-gwt",
    )

    sim.register_ims_package(ims_gwt, [gwt_name])

    dis_gwt = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwt_name}.dis",
    )

    ic = flopy.mf6.ModflowGwtic(
        gwt,
        strt=0,
        filename=f"{gwt_name}.ic",
    )

    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM", filename=f"{gwt_name}.adv")

    porosity = 0.25
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwt_name}.mst")

    ssm = flopy.mf6.ModflowGwtssm(
        gwt,
        filename=f"{gwt_name}.ssm",
    )

    oc_gwt = flopy.mf6.ModflowGwtoc(
        gwt,
        printrecord=[("budget", "all")],
        filename=f"{gwt_name}.oc",
    )

    sftpackagedata = []
    for irno in range(nreaches):
        t = (irno, 0.0)
        sftpackagedata.append(t)

    sft_spd = [(0, "INFLOW", INFLOW_CONC)]
    sft_spd += [(i, "RAINFALL", 0.0) for i in range(ncol)]

    sft_budgetcsv_filerecord = f"{gwt_name}.sft.budget.csv"

    sft = flopy.mf6.ModflowGwtsft(
        gwt,
        boundnames=True,
        budgetcsv_filerecord=sft_budgetcsv_filerecord,
        packagedata=sftpackagedata,
        reachperioddata=sft_spd,
        pname="SFR-1",
        filename=f"{gwt_name}.sft",
    )
    sft_obs = [(f"CONC{n + 1:02d}", "CONCENTRATION", (n,)) for n in range(nreaches)]
    sft_obs += [(f"STORAGE{n + 1:02d}", "STORAGE", (n,)) for n in range(nreaches)]
    sft_obs = {f"{gwt_name}.sft.obs.csv": sft_obs}
    sft.obs.initialize(continuous=sft_obs)

    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwf_name,
        exgmnameb=gwt_name,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_mass_storage_change(times, volumes, concentrations, answers, verbose=False):
    ntimes = len(volumes)
    times = np.array([0.0] + times.tolist())
    volumes = np.array([0.0] + volumes.tolist())
    concentrations = np.array([0.0] + concentrations.tolist())
    mass_storage_changes = np.zeros(ntimes, dtype=float)
    for idx in range(ntimes):
        dt = times[idx + 1] - times[idx]
        v0, v1 = volumes[idx], volumes[idx + 1]
        if v1 == 0.0:
            v0 = 0.0
        mass_storage_changes[idx] = (
            concentrations[idx] * v0 - concentrations[idx + 1] * v1
        ) / dt
    diff = mass_storage_changes - answers

    if verbose:
        print("Mass storage changes:", mass_storage_changes)
        print("Answers:", answers)
        print("Difference:", diff)

    if not np.allclose(diff, 0.0):
        print(f"Mass storage change check failed with max error {np.abs(diff).max()}")
        success = False
    else:
        print("Mass storage change check passed.")
        success = True
    return success


def check_output(idx, test):
    ws = test.workspace
    name = test.name
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws, sim_name=name)
    gwtname = f"{name}_gwt"
    gwt = sim.get_model(gwtname)

    sft_obs_data = gwt.sft.output.obs().get_data()
    for n in range(1, 3):
        tag1 = f"V{n:02d}"
        tag2 = f"CONC{n:02d}"
        tag3 = f"STORAGE{n:02d}"
        success = check_mass_storage_change(
            sft_obs_data["totim"],
            SFR_VOLUME_TS[tag1],
            sft_obs_data[tag2],
            sft_obs_data[tag3],
            verbose=True,
        )

        assert success, f"mass storage change check failed for reach {n}"


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

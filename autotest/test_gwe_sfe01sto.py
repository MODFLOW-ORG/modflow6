"""
Simple one-cell model with two SFR reaches connected in series. Purpose is to
test that the storage term has been successfully corrected when (energy)
transport is active. SFR inflows and rainfall in the SFR network were specified
to result in a defined stage. Initially, the inflow to the reaches is only
from the upstream boundary. However, later in the simulation rainfall becomes
the primary source of inflow to each reach. There is no flow between the stream
and the aquifer.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["sfe_01sto"]

rhk = 0.0
mannings = 0.01
rlen = 1.0
rwid = 1.0
rarea = rlen * rwid
slope = 0.001
inflow_temp = 1.0

sfr_depth = np.array(
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
inflow_depth = np.array(
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

inflow_total = (1.0 / mannings) * rwid * sfr_depth ** (5 / 3) * np.sqrt(slope)
inflow_rate = (1.0 / mannings) * rwid * inflow_depth ** (5 / 3) * np.sqrt(slope)
recharge_rate = inflow_total - inflow_rate
recharge_depth = recharge_rate / rarea
q_mp = inflow_rate + 0.5 * recharge_rate
sfr_depth_mp = (q_mp * mannings / (rwid * np.sqrt(slope))) ** (3 / 5)

dtype = [("V01", float), ("V02", float)]
sfr_volume_ts = np.array(
    [(d0 * rarea, d1 * rarea) for d0, d1 in zip(sfr_depth_mp, sfr_depth)], dtype=dtype
)

# Set some static heat transport related model parameter values
cpw = 4183.0
rhow = 1000.0
cps = 760.0
rhos = 1500.0

# For SFE
K_therm_strmbed = 0.0
rbthcnd = 0.0001


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]

    nlay, nrow, ncol = 1, 1, 1
    delr = rlen
    delc = 1.0
    top = 0.0
    botm = [-100.0]
    nper = len(inflow_total)
    perioddata = [(1.0, 1, 1.0)] * nper

    # <ifno> <cellid> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> <man> <ncon> <ustrf> <ndv>
    conn_data = [[0, -1], [1, 0]]
    package_data = [
        (0, (-1, -1, -1), rlen, rwid, slope, top, 1.0, rhk, mannings, 1, 1.0, 0),
        (1, (-1, -1, -1), rlen, rwid, slope, top, 1.0, rhk, mannings, 1, 1.0, 0),
    ]
    nreaches = len(package_data)

    sfr_dict = {}
    for n in range(nper):
        sfr_spd = [
            (0, "inflow", inflow_rate[n]),
            (0, "rainfall", recharge_depth[n]),
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

    # gwe model
    gwe_name = f"{name}_gwe"
    gwe = flopy.mf6.ModflowGwe(
        sim, modelname=gwe_name, model_nam_file=f"{gwe_name}.nam", save_flows=True
    )

    ims_gwe = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
        linear_acceleration="bicgstab",
        outer_maximum=200,
        inner_maximum=100,
        filename=f"{gwe_name}.ims",
        pname="ims-gwe",
    )

    sim.register_ims_package(ims_gwe, [gwe_name])

    dis_gwe = flopy.mf6.ModflowGwedis(
        gwe,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        top=top,
        botm=botm,
        idomain=1,
        filename=f"{gwe_name}.dis",
    )

    ic = flopy.mf6.ModflowGweic(
        gwe,
        strt=0,
        filename=f"{gwe_name}.ic",
    )

    adv = flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM", filename=f"{gwe_name}.adv")

    porosity = 0.25
    est = flopy.mf6.ModflowGweest(
        gwe,
        heat_capacity_water=cpw,
        density_water=rhow,
        porosity=porosity,
        heat_capacity_solid=cps,
        density_solid=rhos,
        filename=f"{gwe_name}.est",
    )

    ssm = flopy.mf6.ModflowGwessm(
        gwe,
        filename=f"{gwe_name}.ssm",
    )

    oc_gwe = flopy.mf6.ModflowGweoc(
        gwe,
        printrecord=[("budget", "all")],
        filename=f"{gwe_name}.oc",
    )

    sfepackagedata = []
    for irno in range(nreaches):
        t = (irno, 0.0, K_therm_strmbed, rbthcnd)
        sfepackagedata.append(t)

    sfe_spd = [(0, "INFLOW", inflow_temp)]
    sfe_spd += [(i, "RAINFALL", 0.0) for i in range(ncol)]

    sfe_budgetcsv_filerecord = f"{gwe_name}.sfe.budget.csv"

    sfe = flopy.mf6.ModflowGwesfe(
        gwe,
        boundnames=True,
        budgetcsv_filerecord=sfe_budgetcsv_filerecord,
        packagedata=sfepackagedata,
        reachperioddata=sfe_spd,
        pname="SFR-1",
        filename=f"{gwe_name}.sfe",
    )
    sfe_obs = [(f"TEMP{n + 1:02d}", "TEMPERATURE", (n,)) for n in range(nreaches)]
    sfe_obs += [(f"STORAGE{n + 1:02d}", "STORAGE", (n,)) for n in range(nreaches)]
    sfe_obs = {f"{gwe_name}.sfe.obs.csv": sfe_obs}
    sfe.obs.initialize(continuous=sfe_obs)

    gwfgwe = flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwf_name,
        exgmnameb=gwe_name,
        filename=f"{name}.gwfgwe",
    )

    return sim, None


def check_energy_storage_change(times, volumes, temperatures, answers, verbose=False):
    ntimes = len(volumes)
    times = np.array([0.0] + times.tolist())
    volumes = np.array([0.0] + volumes.tolist())
    temperatures = np.array([0.0] + temperatures.tolist())
    energy_storage_changes = np.zeros(ntimes, dtype=float)
    for idx in range(ntimes):
        dt = times[idx + 1] - times[idx]
        v0, v1 = volumes[idx], volumes[idx + 1]
        if v1 == 0.0:
            v0 = 0.0
        energy_storage_changes[idx] = (
            temperatures[idx] * v0 * cpw * rhow
            - temperatures[idx + 1] * v1 * cpw * rhow
        ) / dt
    diff = energy_storage_changes - answers

    if verbose:
        print("Energy storage changes:", energy_storage_changes)
        print("Answers:", answers)
        print("Difference:", diff)

    if not np.allclose(diff, 0.0):
        print(f"Energy storage change check failed with max error {np.abs(diff).max()}")
        success = False
    else:
        print("Energy storage change check passed.")
        success = True
    return success


def check_output(idx, test):
    ws = test.workspace
    name = test.name
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws, sim_name=name)
    gwename = f"{name}_gwe"
    gwe = sim.get_model(gwename)

    sfe_obs_data = gwe.sfe.output.obs().get_data()
    for n in range(1, 3):
        tag1 = f"V{n:02d}"
        tag2 = f"TEMP{n:02d}"
        tag3 = f"STORAGE{n:02d}"
        success = check_energy_storage_change(
            sfe_obs_data["totim"],
            sfr_volume_ts[tag1],
            sfe_obs_data[tag2],
            sfe_obs_data[tag3],
            verbose=True,
        )

        assert success, f"Energy storage change check failed for reach {n}"


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

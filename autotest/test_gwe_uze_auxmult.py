"""
A simple test of confirm that when the auxmultname option is active in UZF,
UZE will stop with an error msg.  UZF objects must have the same area as the
host cell.  If this condition is violated, the code exits with an
appropriate message. Two variations of the same test problem are scripted; the
first test combines GWF and GWE into a single simulation and should error out.
The second test splits GWF and GWE into two different simulations connected
via FMI and the simulation containing the GWE model also should error out.
"""

import flopy
import numpy as np

cases = ["uzauxmlt-err", "uzauxmlt-fmi"]

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97

uzarea_data = {
    "uzf_pkdat": [
        [0, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
        [1, (0, 0, 3), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf01"],
        [2, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
        [3, (0, 0, 4), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf02"],
        [4, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
        [5, (0, 0, 5), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf03"],
        [6, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
        [7, (0, 0, 6), 1, -1, 1.0, 1e-5, 0.2, 0.3, 0.25, 3.5, "uzf04"],
    ],
    "auxmultval": 0.5,
}

nlay, nrow, ncol = 3, 1, 10
nper = 1
perlen = [1.0]
nstp = [1]
tsmult = [1.0]

delr = 100.0
delc = 100.0
strt = -9.0

strt_temp = 1.0
dispersivity = 1.0
cpw = 4180.0
rhow = 1000.0
cps = 800.0
rhos = 2500.0
prsity = 0.2

idomain = np.ones((nlay, nrow, ncol))
idomain[0, 0, 0] = 0
idomain[0, 0, ncol - 1] = 0


def build_gwf(gwf, gwfname):
    flopy.mf6.ModflowGwfdis(
        gwf,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=[-10.0, -20.0, -30.0],
        idomain=idomain,
    )

    # initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=1,
        k=1.0e-4,
        k22=1.0e-4,
        k33=1.0e-5,
    )

    # aquifer storage
    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=1e-5,
        sy=prsity,
        transient=True,
    )

    # chd files
    chdval = -9.0
    iface = 0
    temperature = 10.0
    chdspd = {
        0: [
            [(2, 0, 0), chdval, iface, temperature, "object0"],
            [(2, 0, ncol - 1), chdval, iface, temperature, "object0"],
        ]
    }
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary=["iface", "temperature"],
        boundnames=True,
        print_input=True,
        save_flows=True,
        stress_period_data=chdspd,
        pname="CHD",
    )

    # transient uzf info
    # ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
    uz_pkdat = uzarea_data["uzf_pkdat"]

    finf = 1.0
    extdp = 0.0
    extwc = 0.0
    pet = 0.0
    zero = 0.0
    auxmultval = uzarea_data["auxmultval"]
    uzf_spd = {
        0: [
            [i, finf, pet, extdp, extwc, zero, zero, zero, auxmultval]
            for i in np.arange(len(uz_pkdat))
        ]
    }

    flopy.mf6.ModflowGwfuzf(
        gwf,
        print_input=True,
        print_flows=True,
        save_flows=True,
        boundnames=True,
        ntrailwaves=7,
        nwavesets=40,
        auxiliary="multiplier",
        auxmultname="multiplier",
        package_convergence_filerecord=f"{gwfname}.UzfConvergence.csv",
        wc_filerecord=f"{gwfname}.wc",
        nuzfcells=len(uz_pkdat),
        packagedata=uz_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{gwfname}.uzf.bud",
        pname="UZF",
        filename=f"{gwfname}.uzf",
    )

    # output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwfname}.oc",
    )

    return gwf


def build_gwe(gwe, gwename, fmi=False):
    # Instantiating MODFLOW 6 transport discretization package
    flopy.mf6.ModflowGwedis(
        gwe,
        nogrb=True,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=[-10.0, -20.0, -30.0],
        pname="DIS",
        idomain=idomain,
        filename=f"{gwename}.dis",
    )

    # Instantiating MODFLOW 6 transport initial concentrations
    flopy.mf6.ModflowGweic(
        gwe,
        strt=strt_temp,
        pname="IC",
        filename=f"{gwename}.ic",
    )

    # Instantiating MODFLOW 6 transport advection package
    flopy.mf6.ModflowGweadv(gwe, scheme="TVD", pname="ADV", filename=f"{gwename}.adv")

    # Instantiating MODFLOW 6 transport dispersion package
    flopy.mf6.ModflowGwecnd(
        gwe,
        xt3d_off=False,
        alh=dispersivity,
        ath1=dispersivity,
        ktw=0.5918 * 86400,
        kts=0.2700 * 86400,
        pname="CND",
        filename=f"{gwename}.cnd",
    )

    # Instantiating MODFLOW 6 transport mass storage package
    flopy.mf6.ModflowGweest(
        gwe,
        save_flows=True,
        porosity=prsity,
        heat_capacity_water=cpw,
        density_water=rhow,
        heat_capacity_solid=cps,
        density_solid=rhos,
        pname="EST",
        filename=f"{gwename}.est",
    )

    # Instantiating MODFLOW 6 transport source-sink mixing package
    srctype = "AUX"
    auxname = "TEMPERATURE"
    pname = ["CHD"]
    # Inpput to SSM is: <pname> <srctype> <auxname>
    sources = [[itm, srctype, auxname] for itm in pname]

    flopy.mf6.ModflowGwessm(
        gwe,
        sources=sources,
        pname="SSM",
        filename=f"{gwename}.ssm",
    )

    # Instantiating MODFLOW 6 energy transport source-sink mixing package
    uz_pkdat = uzarea_data["uzf_pkdat"]
    uzepackagedata = [(ct, 1.0) for ct, iuz in enumerate(uz_pkdat)]
    uzeperioddata = {
        0: [[ct, "INFILTRATION", 1.0] for ct, itm in enumerate(uz_pkdat)],
    }

    flopy.mf6.ModflowGweuze(
        gwe,
        flow_package_name="UZF",
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_temperature=True,
        temperature_filerecord=gwename + ".uze.bin",
        budget_filerecord=gwename + ".uze.bud",
        packagedata=uzepackagedata,
        uzeperioddata=uzeperioddata,
        pname="UZE",
        filename=f"{gwename}.uze",
    )

    # Instantiate MODFLOW 6 heat transport output control package
    flopy.mf6.ModflowGweoc(
        gwe,
        pname="OC",
        budget_filerecord=f"{gwename}.cbc",
        temperature_filerecord=f"{gwename}.ucn",
        temperatureprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("TEMPERATURE", "ALL"), ("BUDGET", "ALL")],
        filename=f"{gwename}.oc",
    )

    if fmi:
        gwfname = gwename.replace("gwe", "gwf")
        pd = [
            ("GWFHEAD", gwfname + ".hds", None),
            ("GWFBUDGET", gwfname + ".cbc", None),
            ("UZF", gwfname + ".uzf.bud", None),
        ]
        fmi = flopy.mf6.ModflowGwefmi(gwe, packagedata=pd)

    return gwe


def build_single_sim(idx, ws, exe):
    name = cases[idx]
    gwfname = "gwf-" + name
    gwename = "gwe-" + name

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)

    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions="NEWTON", save_flows=True
    )

    # build out gwf model
    gwf = build_gwf(gwf, gwfname)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    # ----------------------------------
    # Instantiating MODFLOW 6 GWE model
    # ----------------------------------
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    gwe.name_file.save_flows = True

    # build out gwe model
    gwe = build_gwe(gwe, gwename)

    imsgwe = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    # Instantiate Gwf-Gwe Exchange package
    flopy.mf6.ModflowGwfgwe(
        sim,
        exgtype="GWF6-GWE6",
        exgmnamea=gwfname,
        exgmnameb=gwename,
        filename=f"{gwename}.gwfgwe",
    )

    return sim


def run_single_sim(dir, exe):
    idx = 0
    sim = build_single_sim(idx, dir, exe)
    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)
    errmsg = f"simulation should terminate with error message, but did not.\n{buff}"
    assert not success, errmsg


def build_gwf_sim_only(idx, dir, exe):
    name = cases[idx]
    gwfname = "gwf-" + name

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=dir)

    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=gwfname, newtonoptions="NEWTON", save_flows=True
    )

    # build out gwf model
    gwf = build_gwf(gwf, gwfname)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="MODERATE",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    return sim


def build_gwe_sim_only(idx, dir, exe):
    name = cases[idx]
    gwename = "gwe-" + name

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=dir)

    # create tdis package
    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwe = flopy.mf6.ModflowGwe(
        sim, modelname=gwename, model_nam_file=f"{gwename}.nam", save_flows=True
    )
    gwe = build_gwe(gwe, gwename, fmi=True)

    imsgwe = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(imsgwe, [gwe.name])

    return sim


def run_separate_sims(dir, exe):
    idx = 1
    sim = build_gwf_sim_only(idx, dir, exe)
    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)

    # So long as the flow simulation ran successfully, then run GWE'
    errmsg = "GWF only simulation did not run and should have"
    assert success, errmsg
    sim = build_gwe_sim_only(idx, dir, exe)
    sim.write_simulation()
    success, buff = sim.run_simulation(silent=False)
    errmsg = f"GWE simulation should terminate with error message, but did not.\n{buff}"
    assert not success, errmsg


def test_auxmultname(function_tmpdir, targets):
    # Whether running the GWF & GWE model in a single sim or
    # in separate sims via FMI, both should error out, which is
    # how the assertion statements are constructed.
    # Start with check where both GWF and GWE are in same simulation
    run_single_sim(str(function_tmpdir), targets["mf6"])

    # Next, separate GWF and GWE models into 2 different simulations
    # and ensure MF quits with error message
    run_separate_sims(str(function_tmpdir), targets["mf6"])

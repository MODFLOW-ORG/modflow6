"""
Test for checking maw observation input.  The following observation types:
'maw' and 'conductance,' require that ID2 be provided when
ID is an integer corresponding to a well number and not BOUNDNAME.
See table in MAW Package section of mf6io.pdf for an explanation of ID,
ID2, and Observation Type.
"""

import os

import flopy

newtonoptions = [None, "NEWTON", "NEWTON UNDER_RELAXATION"]
cases = "maw_obs"


def build_model(dir, exe):
    nlay, nrow, ncol = 1, 1, 3
    nper = 3
    perlen = [1.0, 1.0, 1.0]
    nstp = [1, 1, 1]
    tsmult = [1.0, 1.0, 1.0]
    lenx = 300.0
    delr = delc = lenx / float(nrow)
    strt = 100.0
    hnoflo = 1e30
    hdry = -1e30
    hk = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 1.0
    krylov = ["CG", "BICGSTAB", "BICGSTAB"]

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
    )
    gwf.name_file.newtonoptions = newtonoptions[0]

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=krylov[0],
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )
    sim.register_ims_package(ims, [gwf.name])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=100.0,
        botm=0.0,
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=1,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=0.0,
        sy=0.1,
        steady_state={0: True},
        # transient={1: False},
        filename=f"{name}.sto",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 100.0])
    chdlist0.append([(0, 0, 2), 100.0])

    chdlist1 = []
    chdlist1.append([(0, 0, 0), 25.0])
    chdlist1.append([(0, 0, 2), 25.0])

    chdspdict = {0: chdlist0, 1: chdlist1, 2: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # wel files
    # wel = flopy.mf6.ModflowGwfwel(gwf, print_input=True, print_flows=True,
    #                              maxbound=len(ws),
    #                              periodrecarray=wd6,
    #                              save_flows=False)
    # MAW
    opth = f"{name}.maw.obs"
    wellbottom = 50.0
    wellrecarray = [[0, 0.1, wellbottom, 100.0, "THIEM", 1]]
    wellconnectionsrecarray = [[0, 0, (0, 0, 1), 100.0, wellbottom, 1.0, 0.1]]
    wellperiodrecarray = [[0, "rate", 0.0]]
    mawo_dict = {}
    mawo_dict["maw_obs.csv"] = [("mh1", "head", 1), ("mawgw", "maw", 1)]
    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        filename=f"{name}.maw",
        print_input=True,
        print_head=True,
        print_flows=True,
        save_flows=True,
        observations=mawo_dict,
        packagedata=wellrecarray,
        connectiondata=wellconnectionsrecarray,
        perioddata=wellperiodrecarray,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim


def test_mf6model(function_tmpdir, targets):
    mf6 = targets["mf6"]
    sim = build_model(str(function_tmpdir), mf6)
    sim.write_simulation()
    sim.run_simulation()

    # ensure that the error msg is contained in the mfsim.lst file
    f = open(str(function_tmpdir / "mfsim.lst"), "r")
    lines = f.readlines()
    error_count = 0
    expected_msg = False
    for line in lines:
        if "ID2 (icon) is missing" in line:
            expected_msg = True
            error_count += 1

    assert error_count == 1, (
        "error count = " + str(error_count) + ", but should equal 1"
    )

    # fix the error and attempt to rerun model
    orig_fl = str(function_tmpdir / (cases + ".maw.obs"))
    new_fl = str(function_tmpdir / (cases + ".maw.obs.new"))
    sr = open(orig_fl, "r")
    sw = open(new_fl, "w")

    lines = sr.readlines()
    error_free_line = "  mawgw  maw  1  1\n"
    for line in lines:
        if " maw " in line:
            sw.write(error_free_line)
        else:
            sw.write(line)

    sr.close()
    sw.close()

    # delete original and replace with corrected lab obs input
    os.remove(orig_fl)
    os.rename(new_fl, orig_fl)

    # rerun the model, should be no errors
    success, buff = sim.run_simulation()
    assert success, buff

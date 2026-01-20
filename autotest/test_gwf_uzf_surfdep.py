import flopy

testname = "uzf_3lay_srfdchk"

iuz_cell_dict = {}
cell_iuz_dict = {}


def build_model(dir, exe):
    nlay, nrow, ncol = 3, 1, 10
    nper = 1
    perlen = [20.0]
    nstp = [10]
    tsmult = len(perlen) * [1.0]

    delr = 1.0
    delc = 1.0
    strt = -25
    botm = [
        [-5.0, -4.0, -3.0, -3.0, -2.0, -5.0, -4.0, -3.0, -3.0, -2.0],
        [-20, -20, -20, -20, -20, -20, -20, -20, -20, -20],
        [-30, -30, -30, -30, -30, -30, -30, -30, -30, -30],
    ]

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-9, 1e-3, 0.97

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = testname

    # build MODFLOW 6 files
    ws = dir
    sim = flopy.mf6.MFSimulation(sim_name=name, version="mf6", exe_name=exe, sim_ws=ws)

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions="NEWTON", save_flows=True
    )

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

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=0.0,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=100.0, k33=10)

    # aquifer storage
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient=True)

    # chd files
    chdval = -3.0
    chdspd = {0: [[(2, 0, 0), chdval]]}
    chd = flopy.mf6.ModflowGwfchd(gwf, print_flows=True, stress_period_data=chdspd)

    # transient uzf info
    # ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
    uzf_pkdat = [
        [0, (0, 0, 1), 1, 8, 6, 1, 0.05, 0.35, 0.05, 4, "uzf01"],
        [1, (0, 0, 2), 1, 9, 6, 1, 0.05, 0.35, 0.05, 4, "uzf02"],
        [2, (0, 0, 3), 1, 10, 6, 1, 0.05, 0.35, 0.05, 4, "uzf03"],
        [3, (0, 0, 4), 1, 11, 6, 1, 0.05, 0.35, 0.05, 4, "uzf04"],
        [4, (0, 0, 5), 1, 12, 6, 1, 0.05, 0.35, 0.05, 4, "uzf05"],
        [5, (0, 0, 6), 1, 13, 6, 1, 0.05, 0.35, 0.05, 4, "uzf06"],
        [6, (0, 0, 7), 1, 14, 6, 1, 0.05, 0.35, 0.05, 4, "uzf07"],
        [7, (0, 0, 8), 1, 15, 6, 1, 0.05, 0.35, 0.05, 4, "uzf08"],
        [8, (1, 0, 1), 0, 16, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf08"],
        [9, (1, 0, 2), 0, 17, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf10"],
        [10, (1, 0, 3), 0, 18, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf11"],
        [11, (1, 0, 4), 0, 19, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf12"],
        [12, (1, 0, 5), 0, 20, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf13"],
        [13, (1, 0, 6), 0, 21, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf14"],
        [14, (1, 0, 7), 0, 22, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf15"],
        [15, (1, 0, 8), 0, 23, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf16"],
        [16, (2, 0, 1), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf17"],
        [17, (2, 0, 2), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf18"],
        [18, (2, 0, 3), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf19"],
        [19, (2, 0, 4), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf20"],
        [20, (2, 0, 5), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf21"],
        [21, (2, 0, 6), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf22"],
        [22, (2, 0, 7), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf23"],
        [23, (2, 0, 8), 0, -1, 0.1, 1, 0.05, 0.35, 0.05, 4, "uzf24"],
    ]

    for itm in uzf_pkdat:
        iuz_cell_dict.update({itm[0]: (itm[1][0], itm[1][1], itm[1][2])})
        cell_iuz_dict.update({(itm[1][0], itm[1][1], itm[1][2]): itm[0]})

    extdp = 15.0
    pet = 0.001
    zero = 0.0
    uzf_spd = {
        0: [
            [0, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [1, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [2, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [3, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [4, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [5, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [6, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [7, 0.01, pet, extdp, 7.0e-02, zero, zero, zero],
            [8, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [9, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [10, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [11, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [12, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [13, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [14, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [15, zero, pet, extdp, 7.0e-02, zero, zero, zero],
        ]
    }

    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        print_flows=True,
        save_flows=True,
        simulate_et=True,
        simulate_gwseep=True,
        linear_gwet=True,
        boundnames=True,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{name}.uzf.bud",
        filename=f"{name}.uzf",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        filename=f"{name}.oc",
    )

    return sim


def test_mf6model(function_tmpdir, targets):
    # build and run the test model
    mf6 = targets["mf6"]
    sim = build_model(str(function_tmpdir), mf6)
    sim.write_simulation()
    sim.run_simulation()

    # ensure that the error msg is contained in the mfsim.lst file
    f = open(str(function_tmpdir / "mfsim.lst"), "r")
    lines = f.readlines()
    error_count = 0
    for line in lines:
        if "SURFDEP" and "cannot" in line:
            error_count += 1

    assert error_count == 8, "error count = " + str(error_count) + "but should equal 8"

    print("Finished running surfdep check")

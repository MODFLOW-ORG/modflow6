import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["uzf_3lay"]
name = "model"
iuz_cell_dict = {}
cell_iuz_dict = {}

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97


def build_models(idx, test):
    nlay, nrow, ncol = 3, 1, 10
    nper = 5
    perlen = [20.0, 20.0, 20.0, 500.0, 2000.0]
    nstp = [10, 10, 10, 10, 50]
    tsmult = len(perlen) * [1.0]

    delr = 1.0
    delc = 1.0
    strt = -25

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

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
        botm=[-5.0, -20, -30],
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=1, k=100.0, k33=10)

    # aquifer storage
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=1e-5, sy=0.2, transient=True)

    # chd files
    chdval = -3.0
    chdspd = {0: [[(2, 0, 0), chdval], [(2, 0, ncol - 1), chdval]]}
    chd = flopy.mf6.ModflowGwfchd(gwf, print_flows=True, stress_period_data=chdspd)

    # transient uzf info
    # ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
    uzf_pkdat = [
        [0, (0, 0, 1), 1, 8, 1, 1, 0.05, 0.35, 0.05, 4, "uzf01"],
        [1, (0, 0, 2), 1, 9, 1, 1, 0.05, 0.35, 0.05, 4, "uzf02"],
        [2, (0, 0, 3), 1, 10, 1, 1, 0.05, 0.35, 0.05, 4, "uzf03"],
        [3, (0, 0, 4), 1, 11, 1, 1, 0.05, 0.35, 0.05, 4, "uzf04"],
        [4, (0, 0, 5), 1, 12, 1, 1, 0.05, 0.35, 0.05, 4, "uzf05"],
        [5, (0, 0, 6), 1, 13, 1, 1, 0.05, 0.35, 0.05, 4, "uzf06"],
        [6, (0, 0, 7), 1, 14, 1, 1, 0.05, 0.35, 0.05, 4, "uzf07"],
        [7, (0, 0, 8), 1, 15, 1, 1, 0.05, 0.35, 0.05, 4, "uzf08"],
        [8, (1, 0, 1), 0, 16, 1, 1, 0.05, 0.35, 0.05, 4, "uzf08"],
        [9, (1, 0, 2), 0, 17, 1, 1, 0.05, 0.35, 0.05, 4, "uzf10"],
        [10, (1, 0, 3), 0, 18, 1, 1, 0.05, 0.35, 0.05, 4, "uzf11"],
        [11, (1, 0, 4), 0, 19, 1, 1, 0.05, 0.35, 0.05, 4, "uzf12"],
        [12, (1, 0, 5), 0, 20, 1, 1, 0.05, 0.35, 0.05, 4, "uzf13"],
        [13, (1, 0, 6), 0, 21, 1, 1, 0.05, 0.35, 0.05, 4, "uzf14"],
        [14, (1, 0, 7), 0, 22, 1, 1, 0.05, 0.35, 0.05, 4, "uzf15"],
        [15, (1, 0, 8), 0, 23, 1, 1, 0.05, 0.35, 0.05, 4, "uzf16"],
        [16, (2, 0, 1), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf17"],
        [17, (2, 0, 2), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf18"],
        [18, (2, 0, 3), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf19"],
        [19, (2, 0, 4), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf20"],
        [20, (2, 0, 5), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf21"],
        [21, (2, 0, 6), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf22"],
        [22, (2, 0, 7), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf23"],
        [23, (2, 0, 8), 0, -1, 1, 1, 0.05, 0.35, 0.05, 4, "uzf24"],
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
        ],
        1: [
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
        ],
        2: [
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
        ],
        3: [
            [0, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [1, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [2, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [3, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [4, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [5, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [6, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [7, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [8, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [9, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [10, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [11, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [12, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [13, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [14, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [15, zero, pet, extdp, 7.0e-02, zero, zero, zero],
        ],
        4: [
            [0, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [1, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [2, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [3, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [4, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [5, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [6, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [7, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [8, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [9, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [10, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [11, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [12, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [13, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [14, zero, pet, extdp, 7.0e-02, zero, zero, zero],
            [15, zero, pet, extdp, 7.0e-02, zero, zero, zero],
        ],
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
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


def check_output(idx, test):
    ws = test.workspace
    test = flopy.mf6.MFSimulation.load(sim_ws=ws)

    bpth = ws / f"{name}.cbc"
    bobj = flopy.utils.CellBudgetFile(bpth, precision="double")
    bobj.get_unique_record_names()
    # '          STO-SS'
    # '          STO-SY'
    # '    FLOW-JA-FACE'
    # '             CHD'
    # '       UZF-GWRCH'
    # '         UZF-GWD'
    # '        UZF-GWET'

    gwet = bobj.get_data(text="UZF-GWET")
    gwet = np.array(gwet)

    uzpth = os.path.join(ws, f"{name}.uzf.bud")
    uzobj = flopy.utils.CellBudgetFile(uzpth, precision="double")
    uzobj.get_unique_record_names()
    # '    FLOW-JA-FACE'
    # '             GWF'
    # '    INFILTRATION'
    # '         REJ-INF'
    # '            UZET'
    # '         STORAGE'

    uzet = uzobj.get_data(text="UZET")
    uzet = np.array(uzet)

    # convert ndarray to grid dimensions
    tot_stp = 0
    tinfo = test.tdis.perioddata.get_data()
    for itm in tinfo:
        tot_stp += int(itm[1])

    gwet_arr = np.zeros(
        (
            tot_stp,
            test.model.dis.nlay.get_data(),
            test.model.dis.nrow.get_data(),
            test.model.dis.ncol.get_data(),
        )
    )

    uzet_arr = np.zeros(
        (
            tot_stp,
            test.model.dis.nlay.get_data(),
            test.model.dis.nrow.get_data(),
            test.model.dis.ncol.get_data(),
        )
    )

    for tm, dat_stp in enumerate(gwet):
        for x, itm in enumerate(dat_stp):
            iuzno = itm[1] - 1  # convert to 0-based
            lay = iuz_cell_dict[iuzno][0]
            row = iuz_cell_dict[iuzno][1]
            col = iuz_cell_dict[iuzno][2]

            gwet_arr[tm, lay, row, col] = itm[2]

    for tm, dat_stp in enumerate(uzet):
        for x, itm in enumerate(dat_stp):
            iuzno = itm[1] - 1  # convert to 0-based
            lay = iuz_cell_dict[iuzno][0]
            row = iuz_cell_dict[iuzno][1]
            col = iuz_cell_dict[iuzno][2]

            uzet_arr[tm, lay, row, col] = itm[2]

    uzf_strsPerDat = test.model.uzf.perioddata.get_data()
    pet = 0
    for tm in range(tot_stp):
        nstps = 0
        for mstp, itm in enumerate(tinfo):
            nstps += itm[1]
            if tm < nstps:
                break

        for i in range(test.model.dis.nrow.get_data()):
            for j in range(test.model.dis.ncol.get_data()):
                if (0, i, j) in cell_iuz_dict:
                    # For this test, pET only specified in the top layer
                    iuz = cell_iuz_dict[(0, i, j)]
                    for m_row in uzf_strsPerDat[mstp]:
                        if m_row[0] == iuz:
                            pet = float(m_row[2])
                            break
                else:
                    continue

                if pet > 0:
                    aet = (gwet_arr[tm, :, i, j] + uzet_arr[tm, :, i, j]).sum()
                    assert abs(aet) <= pet, (
                        "Actual ET exceeds pET at row "
                        + str(i + 1)
                        + " column "
                        + str(j + 1)
                    )


@pytest.mark.parametrize("idx,name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

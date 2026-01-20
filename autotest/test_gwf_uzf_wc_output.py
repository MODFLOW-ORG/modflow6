import os

import flopy
import flopy.utils.binaryfile as bf
import numpy as np
import pytest
from framework import TestFramework

include_NWT = False

cases = ["uzf_3lay_wc_chk"]

iuz_cell_dict = {}
cell_iuz_dict = {}

nlay, nrow, ncol = 3, 1, 10
nper = 6
perlen = [10.0] * 6
nstp = [10] + [1] * 5
tsmult = [1.0] * len(perlen)

delr = 1.0
delc = 1.0
strt = -25

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-9, 1e-3, 0.97

top = 0.0
botm = [-5.0, -20, -30]

icelltype = 1
k = 100.0
k33 = 10
ss = 1e-5
sy = 0.3

ghbelv1 = -25.5
ghbelv2 = -10.5
ghbcond = 1000.0
nwt_ghb_spdat = {
    0: [[2, 0, 0, ghbelv1, ghbcond], [2, 0, ncol - 1, ghbelv1, ghbcond]],
    3: [[2, 0, 0, ghbelv2, ghbcond], [2, 0, ncol - 1, ghbelv2, ghbcond]],
}

ghbspd = {
    0: [[(2, 0, 0), ghbelv1, ghbcond], [(2, 0, ncol - 1), ghbelv1, ghbcond]],
    3: [[(2, 0, 0), ghbelv2, ghbcond], [(2, 0, ncol - 1), ghbelv2, ghbcond]],
}

# ifno  cellid landflg ivertcn surfdp vks thtr thts thti eps [bndnm]
surfdep1 = 1.0
surfdep2 = 0.001
vks = 0.5
eps = 4.0
thti = 0.06
thtr = 0.05
thts = 0.35
iuzfbnd = [[0, 1, 1, 1, 1, 1, 1, 1, 1, 0]]
uzf_pkdat = [
    [0, (0, 0, 1), 1, 8, surfdep1, vks, thtr, thts, thti, eps, "uzf01"],
    [1, (0, 0, 2), 1, 9, surfdep1, vks, thtr, thts, thti, eps, "uzf02"],
    [2, (0, 0, 3), 1, 10, surfdep1, vks, thtr, thts, thti, eps, "uzf03"],
    [3, (0, 0, 4), 1, 11, surfdep1, vks, thtr, thts, thti, eps, "uzf04"],
    [4, (0, 0, 5), 1, 12, surfdep1, vks, thtr, thts, thti, eps, "uzf05"],
    [5, (0, 0, 6), 1, 13, surfdep1, vks, thtr, thts, thti, eps, "uzf06"],
    [6, (0, 0, 7), 1, 14, surfdep1, vks, thtr, thts, thti, eps, "uzf07"],
    [7, (0, 0, 8), 1, 15, surfdep1, vks, thtr, thts, thti, eps, "uzf08"],
    [8, (1, 0, 1), 0, 16, surfdep2, vks, thtr, thts, thti, eps, "uzf09"],
    [9, (1, 0, 2), 0, 17, surfdep2, vks, thtr, thts, thti, eps, "uzf10"],
    [10, (1, 0, 3), 0, 18, surfdep2, vks, thtr, thts, thti, eps, "uzf11"],
    [11, (1, 0, 4), 0, 19, surfdep2, vks, thtr, thts, thti, eps, "uzf12"],
    [12, (1, 0, 5), 0, 20, surfdep2, vks, thtr, thts, thti, eps, "uzf13"],
    [13, (1, 0, 6), 0, 21, surfdep2, vks, thtr, thts, thti, eps, "uzf14"],
    [14, (1, 0, 7), 0, 22, surfdep2, vks, thtr, thts, thti, eps, "uzf15"],
    [15, (1, 0, 8), 0, 23, surfdep2, vks, thtr, thts, thti, eps, "uzf16"],
    [16, (2, 0, 1), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf17"],
    [17, (2, 0, 2), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf18"],
    [18, (2, 0, 3), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf19"],
    [19, (2, 0, 4), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf20"],
    [20, (2, 0, 5), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf21"],
    [21, (2, 0, 6), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf22"],
    [22, (2, 0, 7), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf23"],
    [23, (2, 0, 8), 0, -1, surfdep2, vks, thtr, thts, thti, eps, "uzf24"],
]

for itm in uzf_pkdat:
    iuz_cell_dict.update({itm[0]: (itm[1][0], itm[1][1], itm[1][2])})
    cell_iuz_dict.update({(itm[1][0], itm[1][1], itm[1][2]): itm[0]})

nuztop = 3
iuzfopt = 1
irunflg = 0
ietflg = 1
iuzfcb1 = 10
iuzfcb2 = 0
ntrail2 = 25
nsets2 = 80

extdp = 3.0
extwc = 0.05
pet0 = 0.0
pet1 = 0.001
pet2 = 0.011
finf0 = 0.125
finf1 = 0.001
finf2 = 0.05
finf3 = 0.01
zero = 0.0
uzf1_finf = {0: finf0, 1: finf1, 2: finf2, 3: finf3, 4: finf1, 5: finf1}
uzf1_pet = {0: pet0, 1: pet1, 2: pet1, 3: pet1, 4: pet2, 5: pet2}
uzf_spd = {
    0: [
        [0, finf0, pet0, extdp, extwc, zero, zero, zero],
        [1, finf0, pet0, extdp, extwc, zero, zero, zero],
        [2, finf0, pet0, extdp, extwc, zero, zero, zero],
        [3, finf0, pet0, extdp, extwc, zero, zero, zero],
        [4, finf0, pet0, extdp, extwc, zero, zero, zero],
        [5, finf0, pet0, extdp, extwc, zero, zero, zero],
        [6, finf0, pet0, extdp, extwc, zero, zero, zero],
        [7, finf0, pet0, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
    1: [
        [0, finf1, pet1, extdp, extwc, zero, zero, zero],
        [1, finf1, pet1, extdp, extwc, zero, zero, zero],
        [2, finf1, pet1, extdp, extwc, zero, zero, zero],
        [3, finf1, pet1, extdp, extwc, zero, zero, zero],
        [4, finf1, pet1, extdp, extwc, zero, zero, zero],
        [5, finf1, pet1, extdp, extwc, zero, zero, zero],
        [6, finf1, pet1, extdp, extwc, zero, zero, zero],
        [7, finf1, pet1, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
    2: [
        [0, finf2, pet1, extdp, extwc, zero, zero, zero],
        [1, finf2, pet1, extdp, extwc, zero, zero, zero],
        [2, finf2, pet1, extdp, extwc, zero, zero, zero],
        [3, finf2, pet1, extdp, extwc, zero, zero, zero],
        [4, finf2, pet1, extdp, extwc, zero, zero, zero],
        [5, finf2, pet1, extdp, extwc, zero, zero, zero],
        [6, finf2, pet1, extdp, extwc, zero, zero, zero],
        [7, finf2, pet1, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
    3: [
        [0, finf3, pet1, extdp, extwc, zero, zero, zero],
        [1, finf3, pet1, extdp, extwc, zero, zero, zero],
        [2, finf3, pet1, extdp, extwc, zero, zero, zero],
        [3, finf3, pet1, extdp, extwc, zero, zero, zero],
        [4, finf3, pet1, extdp, extwc, zero, zero, zero],
        [5, finf3, pet1, extdp, extwc, zero, zero, zero],
        [6, finf3, pet1, extdp, extwc, zero, zero, zero],
        [7, finf3, pet1, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
    4: [
        [0, finf1, pet2, extdp, extwc, zero, zero, zero],
        [1, finf1, pet2, extdp, extwc, zero, zero, zero],
        [2, finf1, pet2, extdp, extwc, zero, zero, zero],
        [3, finf1, pet2, extdp, extwc, zero, zero, zero],
        [4, finf1, pet2, extdp, extwc, zero, zero, zero],
        [5, finf1, pet2, extdp, extwc, zero, zero, zero],
        [6, finf1, pet2, extdp, extwc, zero, zero, zero],
        [7, finf1, pet2, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
    5: [
        [0, finf1, pet2, extdp, extwc, zero, zero, zero],
        [1, finf1, pet2, extdp, extwc, zero, zero, zero],
        [2, finf1, pet2, extdp, extwc, zero, zero, zero],
        [3, finf1, pet2, extdp, extwc, zero, zero, zero],
        [4, finf1, pet2, extdp, extwc, zero, zero, zero],
        [5, finf1, pet2, extdp, extwc, zero, zero, zero],
        [6, finf1, pet2, extdp, extwc, zero, zero, zero],
        [7, finf1, pet2, extdp, extwc, zero, zero, zero],
        [8, zero, zero, zero, extwc, zero, zero, zero],
        [9, zero, zero, zero, extwc, zero, zero, zero],
        [10, zero, zero, zero, extwc, zero, zero, zero],
        [11, zero, zero, zero, extwc, zero, zero, zero],
        [12, zero, zero, zero, extwc, zero, zero, zero],
        [13, zero, zero, zero, extwc, zero, zero, zero],
        [14, zero, zero, zero, extwc, zero, zero, zero],
        [15, zero, zero, zero, extwc, zero, zero, zero],
    ],
}


def build_mf6_model(idx, ws):
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
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
        top=top,
        botm=botm,
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=True, icelltype=icelltype, k=k, k33=k33
    )

    # aquifer storage
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1, ss=ss, sy=sy, transient=True)

    # ghb files
    ghb = flopy.mf6.ModflowGwfghb(gwf, print_flows=True, stress_period_data=ghbspd)

    # transient uzf info
    uzf_obs = {
        f"{name}.uzfobs": [
            ("uzf01_dpth=0.5", "water-content", "uzf01", 0.5),
            ("uzf01_dpth=1.5", "water-content", "uzf01", 1.5),  # Relies on boundnames
            ("uzf01_dpth=2.5", "water-content", "uzf01", 2.5),
            ("uzf01_dpth=3.5", "water-content", "uzf01", 3.5),
            ("uzf01_dpth=4.49", "water-content", "uzf01", 4.49),
            ("uzf09_dpth=7.5", "water-content", "uzf09", 7.5),
            ("uzf17_dpth=1.0", "water-content", "uzf17", 1.0),
        ]
    }
    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        print_flows=True,
        save_flows=True,
        wc_filerecord=name + ".uzfwc.bin",
        simulate_et=True,
        simulate_gwseep=True,
        linear_gwet=True,
        boundnames=True,
        observations=uzf_obs,
        ntrailwaves=ntrail2,
        nwavesets=nsets2,
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

    return sim


def build_mfnwt_model(idx, ws):
    name = cases[idx]

    # build MODFLOW-NWT files
    ws = os.path.join(ws, "mfnwt")

    # Instantiate the MODFLOW model
    mf = flopy.modflow.Modflow(
        modelname=name,
        model_ws=ws,
        version="mfnwt",
        exe_name="mfnwt",
    )

    # Instantiate discretization package
    # units: itmuni=4 (days), lenuni=2 (m)
    flopy.modflow.ModflowDis(
        mf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        nper=nper,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        perlen=perlen,
        nstp=nstp,
        itmuni=4,
        lenuni=2,
        steady=False,
    )

    # Instantiate basic package
    flopy.modflow.ModflowBas(mf, ibound=1, strt=strt)

    # Instantiate layer property flow package
    flopy.modflow.ModflowUpw(mf, laytyp=icelltype, layvka=0, hk=k, vka=k33)

    # Instantiate solver package
    flopy.modflow.ModflowNwt(mf, headtol=1e-8, fluxtol=1)

    # Instantiate link mass-transport package (for writing cell-by-cell
    # water contents)
    flopy.modflow.ModflowLmt(mf, output_file_format="formatted", package_flows=["UZF"])

    # Instantiate general head boundary package
    ghb = flopy.modflow.ModflowGhb(mf, stress_period_data=nwt_ghb_spdat)

    # Instantiate unsaturated-zone flow package
    flopy.modflow.ModflowUzf1(
        mf,
        nuztop=nuztop,
        iuzfopt=iuzfopt,
        irunflg=irunflg,
        ietflg=ietflg,
        ipakcb=iuzfcb1,
        iuzfcb2=iuzfcb2,
        ntrail2=ntrail2,
        nsets=nsets2,
        surfdep=surfdep1,
        iuzfbnd=iuzfbnd,
        vks=vks,
        eps=eps,
        thts=thts,
        thti=thti,
        thtr=thtr,
        finf=uzf1_finf,
        pet=uzf1_pet,
        extdp=extdp,
        extwc=extwc,
        nwt_11_fmt=True,
        specifythti=True,
        specifythtr=True,
    )

    # Instantiate output control (OC) package
    spd = {
        (0, 0): ["save head", "save budget", "print budget"],
        (1, 0): ["save head", "save budget", "print budget"],
        (2, 0): ["save head", "save budget", "print budget"],
        (3, 0): ["save head", "save budget", "print budget"],
        (4, 0): ["save head", "save budget", "print budget"],
        (5, 0): ["save head", "save budget", "print budget"],
    }
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd)

    return mf


def build_models(idx, test):
    # Start by building the MF6 model
    sim = build_mf6_model(idx, test.workspace)

    # Construct MF-NWT model for comparing water contents
    #   Commented out to avoid NWT dependency, but left behind for
    #   local testing if needed in the future.
    if include_NWT:
        mc = build_mfnwt_model(idx, test.workspace)
    else:
        mc = None
    return sim, mc


def check_output(idx, test):
    ws = test.workspace

    # Get the MF6 heads
    fpth = os.path.join(ws, "uzf_3lay_wc_chk.hds")
    hobj = flopy.utils.HeadFile(fpth, precision="double")
    hds = hobj.get_alldata()

    # Get the MF6 water contents
    wcpth = os.path.join(ws, cases[0] + ".uzfwc.bin")
    mf6_wc_obj = bf.HeadFile(wcpth, text="   water-content")

    ckstpkper_wc = mf6_wc_obj.get_kstpkper()

    mf6_wc = []
    for current_kstpkper in ckstpkper_wc:
        wc_rawdat = mf6_wc_obj.get_data(kstpkper=current_kstpkper)
        wc_tmp = np.zeros((nlay, nrow, ncol))
        for i, itm in enumerate(uzf_pkdat):
            lay = itm[1][0]
            row = itm[1][1]
            col = itm[1][2]
            wc_tmp[lay, row, col] = wc_rawdat[0][0][i]

        mf6_wc.append(wc_tmp)

    mf6_wc = np.array(mf6_wc)

    # Retrieve MF-NWT water contents from formatted linker file
    if include_NWT:
        mfnwt_wc = []
        fpth = os.path.join(ws, "mfnwt", "mt3d_link.ftl")
        with open(fpth, "r") as f:
            for line in f:
                if "WATER CONTENT   ".lower() in line.lower():
                    line = next(f)
                    wc_tmp = []
                    while "          10           1           3" not in line:
                        m_arr = line.strip().split()
                        for i, val in enumerate(m_arr):
                            wc_tmp.append(float(val))

                        line = next(f)

                    wc_tmp = np.array(wc_tmp)
                    wc_tmp = wc_tmp.reshape((nlay, nrow, ncol))
                    mfnwt_wc.append(wc_tmp)

        mfnwt_wc = np.array(mfnwt_wc)
    else:
        # The following values "burned in" to script and originally calculated
        # by NWT.
        check_vals = [
            8.7777637e-02,
            0.1155553,
            0.1433329,
            0.1711106,
            0.1988882,
            0.2266658,
            0.2544435,
            0.2621320,
            0.2621320,
            0.2621320,
            0.1421950,
            0.2180280,
            0.1622016,
            0.1203337,
            8.9378804e-02,
        ]
        mfnwt_wc = np.array(check_vals)

    if include_NWT:
        # Compare MF6 results with NWT
        # layer 1
        assert np.allclose(mf6_wc[:, 0, 0, :], mfnwt_wc[:, 0, 0, :], atol=0.01), (
            "mf6 water contents different than NWT"
        )
        # layer 2 (has a mixed condition later in the simulation
        assert np.allclose(mf6_wc[:, 1, 0, :], mfnwt_wc[:, 1, 0, :], atol=0.01), (
            "mf6 water contents different than NWT"
        )
    else:
        assert np.allclose(mf6_wc[:, 0, 0, 1], mfnwt_wc, atol=0.01), (
            "mf6 water contents different than established solution"
        )

    print("Finished running checks")


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

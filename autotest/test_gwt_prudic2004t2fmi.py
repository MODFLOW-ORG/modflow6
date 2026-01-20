"""Tests to ability to run flow model first followed by transport model"""

import os
from os.path import join

import flopy
import numpy as np
import pytest
from conftest import project_root_path

data_path = project_root_path / "autotest" / "data"
model_path = str(data_path / "prudic2004test2")
testgroup = "prudic2004t2fmi"

nlay = 8
nrow = 36
ncol = 23
delr = 405.665
delc = 403.717
top = 100.0
fname = os.path.join(model_path, "bot1.dat")
bot0 = np.loadtxt(fname)
botm = [bot0] + [bot0 - (15.0 * k) for k in range(1, nlay)]
fname = os.path.join(model_path, "idomain1.dat")
idomain0 = np.loadtxt(fname, dtype=int)
idomain = nlay * [idomain0]


def run_flow_model(dir, exe):
    global idomain
    name = "flow"
    gwfname = name
    wsf = join(dir, testgroup, name)
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=wsf, exe_name=exe)
    tdis_rc = [(1.0, 1, 1.0), (365.25 * 25, 1, 1.0)]
    nper = len(tdis_rc)
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, save_flows=True)

    # ims
    hclose = 0.0001
    rclose = 0.1
    nouter = 1000
    ninner = 100
    relax = 0.99
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="CG",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    idomain = dis.idomain.array

    ic = flopy.mf6.ModflowGwfic(gwf, strt=50.0)

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        save_saturation=True,
        icelltype=[1] + 7 * [0],
        k=250.0,
        k33=125.0,
    )

    sto_on = False
    if sto_on:
        sto = flopy.mf6.ModflowGwfsto(
            gwf,
            save_flows=True,
            iconvert=[1] + 7 * [0],
            ss=1.0e-5,
            sy=0.3,
            steady_state={0: True},
            transient={1: False},
        )

    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.bud",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    rch_on = True
    if rch_on:
        rch = flopy.mf6.ModflowGwfrcha(gwf, recharge={0: 4.79e-3}, pname="RCH-1")

    chdlist = []
    fname = os.path.join(model_path, "chd.dat")
    for line in open(fname, "r").readlines():
        ll = line.strip().split()
        if len(ll) == 4:
            k, i, j, hd = ll
            chdlist.append([(int(k) - 1, int(i) - 1, int(j) - 1), float(hd)])
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdlist, pname="CHD-1")

    rivlist = []
    fname = os.path.join(model_path, "riv.dat")
    for line in open(fname, "r").readlines():
        ll = line.strip().split()
        if len(ll) == 7:
            k, i, j, s, c, rb, bn = ll
            rivlist.append(
                [
                    (int(k) - 1, int(i) - 1, int(j) - 1),
                    float(s),
                    float(c),
                    float(rb),
                    bn,
                ]
            )
    rivra = flopy.mf6.ModflowGwfriv.stress_period_data.empty(
        gwf, maxbound=len(rivlist), boundnames=True
    )[0]
    for i, t in enumerate(rivlist):
        rivra[i] = tuple(t)
    fname = os.path.join(model_path, "sfr-packagedata.dat")
    sfrpd = np.genfromtxt(fname, names=True)
    sfrpackagedata = flopy.mf6.ModflowGwfsfr.packagedata.empty(
        gwf, boundnames=True, maxbound=sfrpd.shape[0]
    )
    for name in sfrpackagedata.dtype.names:
        if name in rivra.dtype.names:
            sfrpackagedata[name] = rivra[name]
    for name in sfrpackagedata.dtype.names:
        if name in sfrpd.dtype.names:
            sfrpackagedata[name] = sfrpd[name]
    sfrpackagedata["boundname"] = rivra["boundname"]
    fname = os.path.join(model_path, "sfr-connectiondata.dat")
    with open(fname) as f:
        lines = f.readlines()
    sfrconnectiondata = []
    for line in lines:
        t = line.split()
        c = []
        for v in t:
            i = int(v)
            c.append(i)
        sfrconnectiondata.append(c)
    sfrperioddata = {0: [[0, "inflow", 86400], [18, "inflow", 8640.0]]}

    sfr_obs = {
        (gwfname + ".sfr.obs.csv",): [
            ("reach1leakage", "SFR", "LONGESTRIVERINTHEWORLD1"),
            ("reach2leakage", "SFR", "LONGESTRIVERINTHEWORLD2"),
            ("reach3leakage", "SFR", "LONGESTRIVERINTHEWORLD3"),
            ("reach4leakage", "SFR", "LONGESTRIVERINTHEWORLD4"),
        ],
    }
    sfr_obs["digits"] = 7
    sfr_obs["print_input"] = True
    sfr_obs["filename"] = gwfname + ".sfr.obs"

    sfr_on = True
    if sfr_on:
        sfr = flopy.mf6.ModflowGwfsfr(
            gwf,
            print_stage=True,
            print_flows=True,
            stage_filerecord=gwfname + ".sfr.bin",
            budget_filerecord=gwfname + ".sfr.bud",
            mover=True,
            pname="SFR-1",
            length_conversion=3.28084,
            time_conversion=86400.0,
            boundnames=True,
            nreaches=len(rivlist),
            packagedata=sfrpackagedata,
            connectiondata=sfrconnectiondata,
            perioddata=sfrperioddata,
            observations=sfr_obs,
        )

    fname = os.path.join(model_path, "lakibd.dat")
    lakibd = np.loadtxt(fname, dtype=int)
    lakeconnectiondata = []
    nlakecon = [0, 0]
    lak_leakance = 1.0
    for i in range(nrow):
        for j in range(ncol):
            if lakibd[i, j] == 0:
                continue
            else:
                ilak = lakibd[i, j] - 1
                # back
                if i > 0:
                    if lakibd[i - 1, j] == 0 and idomain[0, i - 1, j]:
                        h = [
                            ilak,
                            nlakecon[ilak],
                            (0, i - 1, j),
                            "horizontal",
                            lak_leakance,
                            0.0,
                            0.0,
                            delc / 2.0,
                            delr,
                        ]
                        nlakecon[ilak] += 1
                        lakeconnectiondata.append(h)
                # left
                if j > 0:
                    if lakibd[i, j - 1] and idomain[0, i, j - 1] == 0:
                        h = [
                            ilak,
                            nlakecon[ilak],
                            (0, i, j - 1),
                            "horizontal",
                            lak_leakance,
                            0.0,
                            0.0,
                            delr / 2.0,
                            delc,
                        ]
                        nlakecon[ilak] += 1
                        lakeconnectiondata.append(h)
                # right
                if j < ncol - 1:
                    if lakibd[i, j + 1] == 0 and idomain[0, i, j + 1]:
                        h = [
                            ilak,
                            nlakecon[ilak],
                            (0, i, j + 1),
                            "horizontal",
                            lak_leakance,
                            0.0,
                            0.0,
                            delr / 2.0,
                            delc,
                        ]
                        nlakecon[ilak] += 1
                        lakeconnectiondata.append(h)
                # front
                if i < nrow - 1:
                    if lakibd[i + 1, j] == 0 and idomain[0, i + 1, j]:
                        h = [
                            ilak,
                            nlakecon[ilak],
                            (0, i + 1, j),
                            "horizontal",
                            lak_leakance,
                            0.0,
                            0.0,
                            delc / 2.0,
                            delr,
                        ]
                        nlakecon[ilak] += 1
                        lakeconnectiondata.append(h)
                # vertical
                v = [
                    ilak,
                    nlakecon[ilak],
                    (1, i, j),
                    "vertical",
                    lak_leakance,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ]
                nlakecon[ilak] += 1
                lakeconnectiondata.append(v)

    lak_obs = {
        (gwfname + ".lak.obs.csv",): [
            ("lake1stage", "STAGE", "lake1"),
            ("lake2stage", "STAGE", "lake2"),
            ("lake1leakage", "LAK", "lake1"),
            ("lake2leakage", "LAK", "lake2"),
        ],
    }
    sfr_obs["digits"] = 7
    sfr_obs["print_input"] = True
    sfr_obs["filename"] = gwfname + ".sfr.obs"

    i, j = np.where(lakibd > 0)
    idomain[0, i, j] = 0
    gwf.dis.idomain.set_data(idomain[0], layer=0, multiplier=[1])

    lakpackagedata = [
        [0, 44.0, nlakecon[0], "lake1"],
        [1, 35.2, nlakecon[1], "lake2"],
    ]
    # <outletno> <lakein> <lakeout> <couttype> <invert> <width> <rough> <slope>
    outlets = [[0, 0, -1, "MANNING", 44.5, 5.000000, 0.03, 0.2187500e-02]]

    lake_on = True
    if lake_on:
        lak = flopy.mf6.ModflowGwflak(
            gwf,
            time_conversion=86400.000,
            print_stage=True,
            print_flows=True,
            stage_filerecord=gwfname + ".lak.bin",
            budget_filerecord=gwfname + ".lak.bud",
            mover=True,
            pname="LAK-1",
            boundnames=True,
            nlakes=2,
            noutlets=len(outlets),
            outlets=outlets,
            packagedata=lakpackagedata,
            connectiondata=lakeconnectiondata,
            observations=lak_obs,
        )

    mover_on = True
    if mover_on:
        maxmvr, maxpackages = 2, 2
        mvrpack = [["SFR-1"], ["LAK-1"]]
        mvrperioddata = [
            ["SFR-1", 5, "LAK-1", 0, "FACTOR", 1.0],
            ["LAK-1", 0, "SFR-1", 6, "FACTOR", 1.0],
        ]
        mvr = flopy.mf6.ModflowGwfmvr(
            gwf,
            maxmvr=maxmvr,
            print_flows=True,
            budget_filerecord=gwfname + ".mvr.bud",
            maxpackages=maxpackages,
            packages=mvrpack,
            perioddata=mvrperioddata,
        )

    sim.write_simulation()
    sim.run_simulation(silent=False)

    fname = gwfname + ".hds"
    fname = os.path.join(wsf, fname)
    hobj = flopy.utils.HeadFile(fname, precision="double")
    head = hobj.get_data()
    hobj.file.close()

    # check node2 values for lake, they should be the lake number (1 or 2)
    # and not the connection number
    budgwflak = gwf.output.budget().get_data(text="LAK")
    budgwflak = budgwflak[0]
    node2sim = budgwflak["node2"]
    node2expected = 67 * [1] + 32 * [2]
    errmsg = f"node2 sim not equal node2 expected\n{node2sim}\n{node2expected}"
    assert np.array_equal(node2sim, node2expected), errmsg

    if lake_on:
        fname = gwfname + ".lak.bin"
        fname = os.path.join(wsf, fname)
        lkstage = None
        if os.path.isfile(fname):
            lksobj = flopy.utils.HeadFile(fname, precision="double", text="stage")
            lkstage = lksobj.get_data().flatten()
            lksobj.file.close()

    if sfr_on:
        fname = gwfname + ".sfr.bin"
        fname = os.path.join(wsf, fname)
        sfstage = None
        if os.path.isfile(fname):
            bobj = flopy.utils.HeadFile(fname, precision="double", text="stage")
            sfstage = bobj.get_data().flatten()
            bobj.file.close()

    if mover_on:
        fname = gwfname + ".mvr.bud"
        fname = os.path.join(wsf, fname)
        bobj = flopy.utils.CellBudgetFile(fname, precision="double")
        ra = bobj.recordarray
        print(ra)
        print(ra.dtype)
        for idx in range(ra.shape[0]):
            d = bobj.get_data(idx=idx)[0]
            if d.shape[0] > 0:
                p1 = ra[idx]["paknam"].decode().strip()
                p2 = ra[idx]["paknam2"].decode().strip()
                print(
                    ra[idx]["kstp"],
                    ra[idx]["kper"],
                    ra[idx]["paknam"],
                    ra[idx]["paknam2"],
                )
                for node, node2, q in d:
                    print(p1, node, p2, node2, q)


def run_transport_model(dir, exe):
    name = "transport"
    gwtname = name
    wst = join(dir, testgroup, name)
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name=exe,
        sim_ws=wst,
        continue_=False,
    )

    tdis_rc = [(1.0, 1, 1.0), (365.25 * 25, 25, 1.0)]
    nper = len(tdis_rc)
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

    # ims
    hclose = 0.001
    rclose = 0.001
    nouter = 50
    ninner = 20
    relax = 0.97
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="DBD",
        under_relaxation_theta=0.7,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=idomain,
    )
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)
    sto = flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="TVD")
    dsp = flopy.mf6.ModflowGwtdsp(gwt, alh=20.0, ath1=2, atv=0.2)
    sourcerecarray = [()]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    cnclist = [
        [(0, 0, 11), 500.0],
        [(0, 0, 12), 500.0],
        [(0, 0, 13), 500.0],
        [(0, 0, 14), 500.0],
        [(1, 0, 11), 500.0],
        [(1, 0, 12), 500.0],
        [(1, 0, 13), 500.0],
        [(1, 0, 14), 500.0],
    ]
    cnc = flopy.mf6.ModflowGwtcnc(
        gwt,
        maxbound=len(cnclist),
        stress_period_data=cnclist,
        save_flows=False,
        pname="CNC-1",
    )

    lktpackagedata = [
        (0, 0.0, 99.0, 999.0, "mylake1"),
        (1, 0.0, 99.0, 999.0, "mylake2"),
    ]
    lktperioddata = [
        (0, "STATUS", "ACTIVE"),
        (1, "STATUS", "ACTIVE"),
    ]
    lkt_obs = {
        (gwtname + ".lkt.obs.csv",): [
            ("lkt1conc", "CONCENTRATION", 1),
            ("lkt2conc", "CONCENTRATION", 2),
            ("lkt1frommvr", "FROM-MVR", (0,)),
            ("lkt2frommvr", "FROM-MVR", (1,)),
            ("lkt1tomvr", "TO-MVR", (0,)),
            ("lkt1bntomvr", "TO-MVR", "mylake1"),
        ],
    }
    lkt_obs["digits"] = 7
    lkt_obs["print_input"] = True
    lkt_obs["filename"] = gwtname + ".lkt.obs"

    lkt_on = True
    if lkt_on:
        lkt = flopy.mf6.modflow.ModflowGwtlkt(
            gwt,
            boundnames=True,
            save_flows=True,
            print_input=True,
            print_flows=True,
            print_concentration=True,
            concentration_filerecord=gwtname + ".lkt.bin",
            budget_filerecord=gwtname + ".lkt.bud",
            packagedata=lktpackagedata,
            lakeperioddata=lktperioddata,
            observations=lkt_obs,
            pname="LAK-1",
            auxiliary=["aux1", "aux2"],
        )

    nreach = 38
    sftpackagedata = []
    for irno in range(nreach):
        t = (irno, 0.0, 99.0, 999.0, f"myreach{irno + 1}")
        sftpackagedata.append(t)

    sftperioddata = [(0, "STATUS", "ACTIVE"), (0, "CONCENTRATION", 0.0)]

    sft_obs = {
        (gwtname + ".sft.obs.csv",): [
            (f"sft{i + 1}conc", "CONCENTRATION", i + 1) for i in range(nreach)
        ]
    }
    # append additional obs attributes to obs dictionary
    sft_obs["digits"] = 7
    sft_obs["print_input"] = True
    sft_obs["filename"] = gwtname + ".sft.obs"

    sft_on = True
    if sft_on:
        sft = flopy.mf6.modflow.ModflowGwtsft(
            gwt,
            boundnames=True,
            save_flows=True,
            print_input=True,
            print_flows=True,
            print_concentration=True,
            concentration_filerecord=gwtname + ".sft.bin",
            budget_filerecord=gwtname + ".sft.bud",
            packagedata=sftpackagedata,
            reachperioddata=sftperioddata,
            observations=sft_obs,
            pname="SFR-1",
            auxiliary=["aux1", "aux2"],
        )

    pd = [
        ("GWFHEAD", "../flow/flow.hds", None),
        ("GWFBUDGET", "../flow/flow.bud", None),
        ("GWFMOVER", "../flow/flow.mvr.bud", None),
        ("LAK-1", "../flow/flow.lak.bud", None),
        ("SFR-1", "../flow/flow.sfr.bud", None),
    ]
    fmi = flopy.mf6.ModflowGwtfmi(gwt, packagedata=pd)

    # mover transport package
    mvt = flopy.mf6.modflow.ModflowGwtmvt(gwt, print_flows=True)

    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        budgetcsv_filerecord=f"{gwtname}.bud.csv",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[
            ("COLUMNS", ncol, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    sim.write_simulation()
    sim.run_simulation()

    fname = gwtname + ".lkt.bin"
    fname = os.path.join(wst, fname)
    bobj = flopy.utils.HeadFile(fname, precision="double", text="concentration")
    lkaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    fname = gwtname + ".sft.bin"
    fname = os.path.join(wst, fname)
    bobj = flopy.utils.HeadFile(fname, precision="double", text="concentration")
    sfaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    # set atol
    atol = 0.05

    # check simulated concentration in lak 1 and 2 sfr reaches
    res_lak1 = lkaconc[:, 0]
    ans_lak1 = [
        -1.73249951e-19,
        5.97398983e-02,
        4.18358112e-01,
        1.48857598e00,
        3.63202585e00,
        6.92925430e00,
        1.11162776e01,
        1.57328143e01,
        2.03088745e01,
        2.45013060e01,
        2.81200704e01,
        3.11132152e01,
        3.34833369e01,
        3.53028319e01,
        3.66693021e01,
        3.76781530e01,
        3.84188513e01,
        3.89615387e01,
        3.93577458e01,
        3.96464993e01,
        3.98598113e01,
        4.00184878e01,
        4.01377654e01,
        4.02288674e01,
        4.02998291e01,
        4.03563314e01,
    ]
    ans_lak1 = np.array(ans_lak1)
    d = res_lak1 - ans_lak1
    msg = f"{res_lak1}\n{ans_lak1}\n{d}"
    assert np.allclose(res_lak1, ans_lak1, atol=atol), msg

    res_sfr3 = sfaconc[:, 30]
    ans_sfr3 = [
        -7.67944651e-23,
        5.11358249e-03,
        3.76169957e-02,
        1.42055634e-01,
        3.72438193e-01,
        7.74112522e-01,
        1.37336373e00,
        2.18151231e00,
        3.19993561e00,
        4.42853144e00,
        5.85660993e00,
        7.46619448e00,
        9.22646330e00,
        1.11069607e01,
        1.30764504e01,
        1.50977917e01,
        1.71329980e01,
        1.91636634e01,
        2.11530199e01,
        2.30688490e01,
        2.48821059e01,
        2.65691424e01,
        2.81080543e01,
        2.94838325e01,
        3.06909748e01,
        3.17352915e01,
    ]
    ans_sfr3 = np.array(ans_sfr3)
    d = res_sfr3 - ans_sfr3
    msg = f"{res_sfr3}\n{ans_sfr3}\n{d}"
    assert np.allclose(res_sfr3, ans_sfr3, atol=atol), msg

    res_sfr4 = sfaconc[:, 37]
    ans_sfr4 = [
        -2.00171747e-20,
        3.55076535e-02,
        2.49465789e-01,
        8.91299656e-01,
        2.18622372e00,
        4.19920114e00,
        6.79501651e00,
        9.72255743e00,
        1.27208739e01,
        1.55989390e01,
        1.82462345e01,
        2.06258607e01,
        2.27255881e01,
        2.45721928e01,
        2.62061367e01,
        2.76640442e01,
        2.89788596e01,
        3.01814571e01,
        3.12842113e01,
        3.22945541e01,
        3.32174210e01,
        3.40539043e01,
        3.48027700e01,
        3.54636082e01,
        3.60384505e01,
        3.65330352e01,
    ]
    ans_sfr4 = np.array(ans_sfr4)
    d = res_sfr4 - ans_sfr4
    msg = f"{res_sfr4}\n{ans_sfr4}\n{d}"
    assert np.allclose(res_sfr4, ans_sfr4, atol=atol), msg

    # make some checks on lake obs csv file
    fname = gwtname + ".lkt.obs.csv"
    fname = os.path.join(wst, fname)
    try:
        tc = np.genfromtxt(fname, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fname}"'
    errmsg = f"to-mvr boundname and outlet number do not match for {fname}"
    assert np.allclose(tc["LKT1TOMVR"], tc["LKT1BNTOMVR"]), errmsg

    # Compare the budget terms from the list file and the budgetcsvfile
    fname = f"{gwtname}.bud.csv"
    fname = os.path.join(wst, fname)
    csvra = np.genfromtxt(fname, dtype=None, names=True, delimiter=",", deletechars="")

    fname = f"{gwtname}.lst"
    fname = os.path.join(wst, fname)
    lst = flopy.utils.Mf6ListBudget(fname, budgetkey="MASS BUDGET FOR ENTIRE MODEL")
    lstra = lst.get_incremental()

    # list file has additional terms, so need to pluck out the following for
    # direct comparison
    imap = (0, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 9, 16, 18)
    success_all = True
    failed_list = []
    for name1, i in zip(csvra.dtype.names, imap):
        name2 = lstra.dtype.names[i]
        if i == 18:
            # percent difference needs to compare small numbers
            success = np.allclose(csvra[name1], lstra[name2], atol=1.0e-7)
        else:
            success = np.allclose(csvra[name1], lstra[name2], rtol=0.01)
        if not success:
            success_all = False
            failed_list.append(name1)
            print(f"Records do not match for {name1} in position {i}")
            for rate1, rate2 in zip(csvra[name1], lstra[name2]):
                print(rate1, rate2)
    assert success_all, f"Comparisons failed for {failed_list}"


@pytest.mark.slow
def test_prudic2004t2fmi(function_tmpdir, targets):
    run_flow_model(str(function_tmpdir), targets["mf6"])
    run_transport_model(str(function_tmpdir), targets["mf6"])

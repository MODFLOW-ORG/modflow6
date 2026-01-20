"""
Tests ATS on the prudic transport model.  With these ATS settings, the
solver should fail on time step 19 in period 2, and should converge on the
second try with a smaller time step.  This test will not pass if the states
are not restored properly for the advanced transport packages when the
failure occurs.
"""

import os
from os.path import join

import flopy
import numpy as np
import pytest
from conftest import project_root_path

data_path = project_root_path / "autotest" / "data"
model_path = str(data_path / "prudic2004test2")
testgroup = "prudic2004t2fmiats"
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
    hclose = 0.001
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
    tdis = flopy.mf6.ModflowTdis(
        sim,
        time_units="DAYS",
        nper=nper,
        perioddata=tdis_rc,
    )

    if True:
        dt0 = 100
        dtmin = 1.0e-5
        dtmax = 10000.0
        dtadj = 2.0
        dtfailadj = 5.0
        ats_filerecord = name + ".ats"
        atsperiod = [
            (1, dt0, dtmin, dtmax, dtadj, dtfailadj),
        ]
        tdis.ats.initialize(
            maxats=len(atsperiod),
            perioddata=atsperiod,
            filename=ats_filerecord,
        )

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
    times = np.array(times)
    ans_times = [
        1.0,
        101.0,
        201.0,
        301.0,
        501.0,
        701.0,
        901.0,
        1101.0,
        1501.0,
        1901.0,
        2301.0,
        2701.0,
        3101.0,
        3501.0,
        3901.0,
        4301.0,
        4701.0,
        5101.0,
        5501.0,
        5581.0,
        5661.0,
        5821.0,
        5981.0,
        6141.0,
        6301.0,
        6621.0,
        6941.0,
        7581.0,
        8221.0,
        8861.0,
        9132.25,
    ]
    ans_times = np.array(ans_times)
    errmsg = "Expected number of total timesteps is different."
    assert times.shape == ans_times.shape, errmsg
    errmsg = f"Times {times} not equal expected times {ans_times}"
    assert np.allclose(times, ans_times)

    # set atol
    atol = 0.05

    # check simulated concentration in lak 1 and 2 sfr reaches
    res_lak1 = lkaconc[:, 0]
    ans_lak1 = [
        -1.7334085635551077e-19,
        -3.187033329925361e-07,
        -1.9287290216857604e-06,
        5.808788660373555e-07,
        0.005591936631026452,
        0.04542773591098022,
        0.1928635682644908,
        0.5690001383534176,
        2.999420704893868,
        7.110019025850782,
        12.268025985205634,
        17.67417138740906,
        22.69808352286938,
        27.00477920391491,
        30.50733530522461,
        33.222437858798955,
        35.25052779893794,
        36.73069024938392,
        37.792799707882,
        37.98952686059535,
        38.171619378463866,
        38.48532541433273,
        38.755615320864834,
        38.98852685483134,
        39.189072004020026,
        39.491640226795035,
        39.71996654913013,
        40.00486884056597,
        40.18758842234358,
        40.309629842366334,
        40.35288988875558,
    ]
    ans_lak1 = np.array(ans_lak1)
    d = res_lak1 - ans_lak1
    msg = f"{res_lak1}\n{ans_lak1}\n{d}"
    assert np.allclose(res_lak1, ans_lak1, atol=atol), msg

    res_sfr3 = sfaconc[:, 30]
    ans_sfr3 = [
        -7.607756096700458e-23,
        -1.3669399889004086e-08,
        -9.199259301584774e-08,
        5.257821481671474e-08,
        0.00039938114816238295,
        0.003355197965954286,
        0.014744417223049597,
        0.04510445881222458,
        0.27877628044373737,
        0.7458007019884897,
        1.4665631236737788,
        2.444128940946191,
        3.6753371432162,
        5.158039470416099,
        6.873092531310018,
        8.780865680435873,
        10.839098670662382,
        13.005360141577922,
        15.230011242915861,
        15.676842775991494,
        16.124955033719825,
        17.022019083397183,
        17.92055413970275,
        18.81985349843973,
        19.717233813700727,
        21.475928749632136,
        23.177014613583257,
        26.206656959204977,
        28.767131611820425,
        30.825804240468084,
        31.611737865014057,
    ]
    ans_sfr3 = np.array(ans_sfr3)
    d = res_sfr3 - ans_sfr3
    msg = f"{res_sfr3}\n{ans_sfr3}\n{d}"
    assert np.allclose(res_sfr3, ans_sfr3, atol=atol), msg

    res_sfr4 = sfaconc[:, 37]
    ans_sfr4 = [
        -2.0041563248238944e-20,
        -1.319932823356964e-07,
        -8.732574723159916e-07,
        7.946044660303284e-07,
        0.00328713663771796,
        0.026759005065411397,
        0.11383700188916444,
        0.33657660610442625,
        1.7925681982241561,
        4.287068052572867,
        7.4770472357597395,
        10.91908508386622,
        14.260965416010272,
        17.315284886138496,
        20.02515736345197,
        22.382282530045032,
        24.423804946506543,
        26.208303689647263,
        27.786040077700754,
        28.093696017884817,
        28.393775338962325,
        28.966213039403637,
        29.515140249737588,
        30.043603904714946,
        30.553203510121918,
        31.501577161886416,
        32.38308331851197,
        33.88529441352757,
        35.12256241115968,
        36.10351180222542,
        36.47615223874849,
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

    # check simulation list file for ats information
    fname = os.path.join(wst, "mfsim.lst")
    with open(fname, "r") as f:
        lines = f.readlines()

    txtlist = [
        (
            "Failed solution for step 19 and period 2 will be retried using "
            "time step of    80.00000"
        ),
        "ATS IS OVERRIDING TIME STEPPING FOR THIS PERIOD",
        "INITIAL TIME STEP SIZE                 (DT0) =    100.0000",
        "MINIMUM TIME STEP SIZE               (DTMIN) =   0.1000000E-04",
        "MAXIMUM TIME STEP SIZE               (DTMAX) =    10000.00",
        "MULTIPLIER/DIVIDER FOR TIME STEP     (DTADJ) =    2.000000",
        "DIVIDER FOR FAILED TIME STEP     (DTFAILADJ) =    5.000000",
    ]
    all_found = True
    for stxt in txtlist:
        msg = f"Checking for string in mfsim.lst: {stxt}"
        found = False
        for line in lines:
            if stxt in line:
                found = True
                break
        if not found:
            msg += " -- NOT FOUND!"
            all_found = False
            print(f"text not found in mfsim.lst: {stxt}")
        print(msg)
    assert all_found, "One or more required text strings not found in mfsim.lst"


@pytest.mark.slow
def test_prudic2004t2fmiats(function_tmpdir, targets):
    run_flow_model(dir=str(function_tmpdir), exe=targets["mf6"])
    run_transport_model(dir=str(function_tmpdir), exe=targets["mf6"])

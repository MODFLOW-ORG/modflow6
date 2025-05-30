"""
Second problem described by Prudic et al (2004)
This problem involves transport through an aquifers, lakes and streams.
It requires the use of the Water Mover Package to send water from a stream,
into a lake, and then back into another stream. Solute is also transport
through the system.
"""

import os

import flopy
import numpy as np
import pytest
from conftest import project_root_path
from framework import TestFramework

cases = ["prudic2004t2"]
data_path = project_root_path / "autotest" / "data"
model_path = data_path / "prudic2004test2"
fname = str(model_path / "lakibd.dat")
lakibd = np.loadtxt(fname, dtype=int)


def build_models(idx, test):
    ws = test.workspace
    name = cases[idx]
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        continue_=False,
        memory_print_option="ALL",
    )

    # number of time steps for period 2 are reduced from 12 * 25 to 25 in
    # order to speed up this autotest
    tdis_rc = [(1.0, 1, 1.0), (365.25 * 25, 25, 1.0)]
    nper = len(tdis_rc)
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname)

    # ims
    nouter, ninner = 1000, 50
    hclose, rclose, relax = 1e-6, 1e-6, 0.99
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_maximum=nouter,
        inner_maximum=ninner,
        outer_dvclose=hclose,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )

    nlay = 8
    nrow = 36
    ncol = 23
    delr = 405.665
    delc = 403.717
    top = 100.0
    fname = str(model_path / "bot1.dat")
    bot0 = np.loadtxt(fname)
    botm = [bot0] + [bot0 - (15.0 * k) for k in range(1, nlay)]
    fname = str(model_path / "idomain1.dat")
    idomain0 = np.loadtxt(fname, dtype=int)
    idomain = nlay * [idomain0]
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
    fname = str(model_path / "chd.dat")
    for line in open(fname, "r").readlines():
        ll = line.strip().split()
        if len(ll) == 4:
            k, i, j, hd = ll
            chdlist.append([(int(k) - 1, int(i) - 1, int(j) - 1), float(hd)])
    chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdlist, pname="CHD-1")

    rivlist = []
    fname = str(model_path / "riv.dat")
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
    sfrpd = np.genfromtxt(model_path / "sfr-packagedata.dat", names=True)
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
    with open(model_path / "sfr-connectiondata.dat") as f:
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
        )

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
            maxpackages=maxpackages,
            packages=mvrpack,
            perioddata=mvrperioddata,
        )

    transport = True
    if transport:
        gwt = flopy.mf6.ModflowGwt(sim, modelname=gwtname)

        # ims
        imsgwt = flopy.mf6.ModflowIms(
            sim,
            print_option="ALL",
            outer_maximum=nouter,
            inner_maximum=ninner,
            outer_dvclose=hclose,
            inner_dvclose=hclose,
            rcloserecord=rclose,
            under_relaxation="DBD",
            under_relaxation_theta=0.7,
            linear_acceleration="bicgstab",
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
        mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.3)
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

        lkt_obs = {}
        for obstype in [
            "CONCENTRATION",
            "STORAGE",
            "CONSTANT",
            "FROM-MVR",
            "TO-MVR",
            "RAINFALL",
            "EVAPORATION",
            "RUNOFF",
            "EXT-INFLOW",
            "WITHDRAWAL",
            "EXT-OUTFLOW",
        ]:
            fname = f"{gwtname}.lkt.obs.{obstype.lower()}.csv"
            obs1 = []
            ncv = 2
            if obstype == "TO-MVR":
                ncv = 1
            obs1 = [(f"lkt{i + 1}", obstype, i + 1) for i in range(ncv)]
            obs2 = [(f"blkt{i + 1}", obstype, f"mylake{i + 1}") for i in range(ncv)]
            lkt_obs[fname] = obs1 + obs2

            # add LKT specific obs
            obstype = "LKT"
            ncv = 2
            nconn = [67, 32]
            fname = f"{gwtname}.lkt.obs.{obstype.lower()}.csv"
            obs1 = [
                (f"lkt{1}-{iconn + 1}", obstype, 1, iconn + 1)
                for iconn in range(nconn[0])  # lake 1
            ] + [
                (f"lkt{2}-{iconn + 1}", obstype, 2, iconn + 1)
                for iconn in range(nconn[1])  # lake 2
            ]
            obs2 = [(f"blkt{i + 1}", obstype, f"mylake{i + 1}") for i in range(ncv)]
            lkt_obs[fname] = obs1 + obs2

        lkt_obs["digits"] = 15
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

        sftpackagedata = []
        for irno in range(sfrpd.shape[0]):
            t = (irno, 0.0, 99.0, 999.0, f"myreach{irno + 1}")
            sftpackagedata.append(t)

        sftperioddata = [(0, "STATUS", "ACTIVE"), (0, "CONCENTRATION", 0.0)]

        sft_obs = {}
        for obstype in [
            "CONCENTRATION",
            "STORAGE",
            "CONSTANT",
            "FROM-MVR",
            "TO-MVR",
            "SFT",
            "RAINFALL",
            "EVAPORATION",
            "RUNOFF",
            "EXT-INFLOW",
            "EXT-OUTFLOW",
        ]:
            fname = f"{gwtname}.sft.obs.{obstype.lower()}.csv"
            obs1 = [(f"sft{i + 1}", obstype, i + 1) for i in range(sfrpd.shape[0])]
            obs2 = [
                (f"bsft{i + 1}", obstype, f"myreach{i + 1}")
                for i in range(sfrpd.shape[0])
            ]
            sft_obs[fname] = obs1 + obs2

        obstype = "FLOW-JA-FACE"
        fname = f"{gwtname}.sft.obs.{obstype.lower()}.csv"
        obs1 = []
        for id1, reach_connections in enumerate(sfrconnectiondata):
            for id2 in reach_connections:
                id2 = abs(id2)
                if id1 != id2:
                    obs1.append((f"sft{id1 + 1}x{id2 + 1}", obstype, id1 + 1, id2 + 1))
        obs2 = [
            (f"bsft{i + 1}", obstype, f"myreach{i + 1}") for i in range(sfrpd.shape[0])
        ]
        sft_obs[fname] = obs1 + obs2

        # append additional obs attributes to obs dictionary
        sft_obs["digits"] = 15
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

        # mover transport package
        mvt = flopy.mf6.modflow.ModflowGwtmvt(gwt)

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

        # GWF GWT exchange
        gwfgwt = flopy.mf6.ModflowGwfgwt(
            sim, exgtype="GWF6-GWT6", exgmnamea=gwfname, exgmnameb=gwtname
        )

    return sim, None


def make_concentration_vs_time(sim):
    print("making plot of concentration versus time...")
    name = sim.name
    ws = sim.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    bobj = gwt.lkt.output.concentration()
    lkaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    bobj = gwt.sft.output.concentration()
    sfaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    times = np.array(times) / 365.0
    plt.plot(times, lkaconc[:, 0], "b-", label="Lake 1")
    plt.plot(times, sfaconc[:, 30], "r-", label="Stream segment 3")
    plt.plot(times, sfaconc[:, 37], "g-", label="Stream segment 4")
    plt.legend()
    plt.ylim(0, 50)
    plt.xlim(0, 25)
    plt.xlabel("TIME, IN YEARS")
    plt.ylabel("SIMULATED BORON CONCENTRATION,\nIN MICROGRAMS PER LITER")
    plt.draw()
    fname = os.path.join(ws, "fig-concentration_vs_time.png")
    print(f"Creating {fname}")
    plt.savefig(fname)

    return


def make_concentration_map(sim):
    print("making concentration map...")

    import matplotlib.pyplot as plt

    levels = [1, 10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]

    name = sim.name
    ws = sim.workspace
    simfp = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = simfp.get_model(gwfname)
    gwt = simfp.get_model(gwtname)
    conc = gwt.output.concentration().get_data()
    lakconc = gwt.lak.output.concentration().get_data().flatten()

    il, jl = np.where(lakibd > 0)
    for i, j in zip(il, jl):
        ilak = lakibd[i, j] - 1
        lake_conc = lakconc[ilak]
        conc[0, i, j] = lake_conc

    fig, axs = plt.subplots(2, 2, figsize=(5, 7), dpi=300, tight_layout=True)

    # plot layers 1, 3, 5, and 8
    for iplot, ilay in enumerate([0, 2, 4, 7]):
        ax = axs.flatten()[iplot]
        pmv = flopy.plot.PlotMapView(model=gwt, ax=ax, layer=ilay)
        # pmv.plot_grid()
        pmv.plot_array(lakibd, masked_values=[0], alpha=0.2)
        pmv.plot_inactive(color_noflow="gray", alpha=0.25)
        # pmv.plot_bc(name="CHD-1", color="blue")
        cs = pmv.contour_array(conc, levels=levels, masked_values=[1.0e30])
        ax.clabel(cs, cs.levels[::1], fmt="%1.0f", colors="b")
        ax.set_title(f"Model layer {ilay + 1}")

    fname = os.path.join(ws, "fig-concentration.png")
    print(f"Creating {fname}")
    plt.savefig(fname)

    return


def check_obs(sim):
    print("checking obs...")
    name = sim.name
    ws = sim.workspace
    sim = flopy.mf6.MFSimulation.load(sim_ws=ws)
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)

    # ensure SFT obs are the same whether specified by
    # boundname or by reach
    csvfiles = gwt.sft.obs.output.obs_names
    for csvfile in csvfiles:
        if ".flow-ja-face.csv" in csvfile:
            continue
        print(f"Checking csv file: {csvfile}")
        conc_ra = gwt.sft.obs.output.obs(f=csvfile).data
        # save a couple entries for comparison with lake
        if ".to-mvr." in csvfile:
            sft6tomvr = conc_ra["BSFT6"]
        if ".from-mvr." in csvfile:
            sft7tomvr = conc_ra["BSFT7"]
        success = True
        for ireach in range(38):
            # print(f"  Checking reach {ireach + 1}")
            is_same = np.allclose(
                conc_ra[f"SFT{ireach + 1}"], conc_ra[f"BSFT{ireach + 1}"]
            )
            if not is_same:
                success = False
                for t, x, y in zip(
                    conc_ra["totim"],
                    conc_ra[f"SFT{ireach + 1}"],
                    conc_ra[f"BSFT{ireach + 1}"],
                ):
                    print(t, x, y)

    # process the sft values and make sure the individual connection rates
    # add up to the boundname rate
    csvfile = "gwt_prudic2004t2.sft.obs.flow-ja-face.csv"
    print(f"Checking csv file: {csvfile}")
    conc_ra = gwt.sft.obs.output.obs(f=csvfile).data
    ntimes = conc_ra.shape[0]
    for ireach in range(38):
        connection_sum = np.zeros(ntimes)
        for column_name in conc_ra.dtype.names:
            if f"SFT{ireach + 1}X" in column_name:
                connection_sum += conc_ra[column_name]
        is_same = np.allclose(connection_sum, conc_ra[f"BSFT{ireach + 1}"])
        if not is_same:
            success = False
            diff = connection_sum - conc_ra[f"BSFT{ireach + 1}"]
            print(
                f"Problem with SFT {ireach + 1}; "
                f"mindiff {diff.min()} and maxdiff {diff.max()}"
            )
            # for itime, (cs, bsft) in enumerate(
            #     zip(connection_sum, conc_ra[f"BSFT{ireach + 1}"])
            # ):
            #     print(itime, cs, bsft)

    assert success, "One or more SFT obs checks did not pass"

    # ensure LKT obs are the same whether specified by
    # boundname or by reach
    csvfiles = gwt.lkt.obs.output.obs_names
    for csvfile in csvfiles:
        if ".lkt.csv" in csvfile:
            continue
        print(f"Checking csv file: {csvfile}")
        conc_ra = gwt.lkt.obs.output.obs(f=csvfile).data
        if ".from-mvr." in csvfile:
            lkt1frommvr = conc_ra["BLKT1"]
        if ".to-mvr." in csvfile:
            lkt1tomvr = conc_ra["BLKT1"]
        success = True
        if ".to-mvr." in csvfile:
            numvalues = 1  # outlet
        else:
            numvalues = 2  # lakes
        for ilake in range(numvalues):
            # print(f"  Checking lake {ilake + 1}")
            is_same = np.allclose(
                conc_ra[f"LKT{ilake + 1}"], conc_ra[f"BLKT{ilake + 1}"]
            )
            if not is_same:
                success = False
                for t, x, y in zip(
                    conc_ra["totim"],
                    conc_ra[f"LKT{ilake + 1}"],
                    conc_ra[f"BLKT{ilake + 1}"],
                ):
                    print(t, x, y)

    # process the lkt values and make sure the individual connection rates
    # add up to the boundname rate
    csvfile = "gwt_prudic2004t2.lkt.obs.lkt.csv"
    print(f"Checking csv file: {csvfile}")
    conc_ra = gwt.lkt.obs.output.obs(f=csvfile).data
    ntimes = conc_ra.shape[0]
    for ilake in [0, 1]:
        connection_sum = np.zeros(ntimes)
        for column_name in conc_ra.dtype.names:
            if f"LKT{ilake + 1}" in column_name and column_name.startswith("LKT"):
                connection_sum += conc_ra[column_name]
        is_same = np.allclose(connection_sum, conc_ra[f"BLKT{ilake + 1}"])
        if not is_same:
            success = False
            print(f"Problem with Lake {ilake + 1}")
            for itime, (cs, blkt) in enumerate(zip(connection_sum, conc_ra["BLKT1"])):
                print(itime, cs, blkt)

    assert success, "One or more LKT obs checks did not pass"

    # check that SFT6 to-mvr is equal to LKT1 from-mvr
    success = True
    is_same = np.allclose(-sft6tomvr, lkt1frommvr, atol=0.1)
    if not is_same:
        success = False
        print("Problem with sft6tomvr comparison to lkt1frommvr")
        for itime, (a, b) in enumerate(zip(-sft6tomvr, lkt1frommvr)):
            print(itime, a, b)

    is_same = np.allclose(-lkt1tomvr, sft7tomvr)
    if not is_same:
        success = False
        print("Problem with lkt1tomvr comparison to sft7tomvr")
        for itime, (a, b) in enumerate(zip(-lkt1tomvr, sft7tomvr)):
            print(itime, a, b)

    assert success, "One or more SFT-LKT obs checks did not pass"


def plot_output(idx, test):
    make_concentration_vs_time(test)
    make_concentration_map(test)


def check_output(idx, test):
    ws = test.workspace
    name = test.name
    gwtname = "gwt_" + name

    check_obs(test)

    fname = gwtname + ".lkt.bin"
    fname = os.path.join(ws, fname)
    bobj = flopy.utils.HeadFile(fname, precision="double", text="concentration")
    lkaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    fname = gwtname + ".sft.bin"
    fname = os.path.join(ws, fname)
    bobj = flopy.utils.HeadFile(fname, precision="double", text="concentration")
    sfaconc = bobj.get_alldata()[:, 0, 0, :]
    times = bobj.times
    bobj.file.close()

    # set atol
    atol = 0.02

    # check simulated concentration in lak 1 and 2 sfr reaches
    res_lak1 = lkaconc[:, 0]
    ans_lak1 = [
        -1.65182115e-19,
        7.68348266e-02,
        5.03641472e-01,
        1.71779241e00,
        4.07354668e00,
        7.62114796e00,
        1.20857818e01,
        1.69944542e01,
        2.18518461e01,
        2.62846600e01,
        3.00831076e01,
        3.31870436e01,
        3.56191870e01,
        3.74752627e01,
        3.88707570e01,
        3.99138871e01,
        4.06883115e01,
        4.12613525e01,
        4.16851085e01,
        4.19995016e01,
        4.22342082e01,
        4.24102803e01,
        4.25433647e01,
        4.26458347e01,
        4.27267325e01,
        4.27923612e01,
    ]
    ans_lak1 = np.array(ans_lak1)
    d = res_lak1 - ans_lak1
    msg = f"{res_lak1} {ans_lak1} {d}"
    assert np.allclose(res_lak1, ans_lak1, atol=atol), msg

    res_sfr3 = sfaconc[:, 30]
    ans_sfr3 = [
        -1.60384621e-22,
        6.41952680e-03,
        4.43870909e-02,
        1.61248207e-01,
        4.12297263e-01,
        8.43111144e-01,
        1.48209168e00,
        2.34273319e00,
        3.42585405e00,
        4.72610177e00,
        6.23101960e00,
        7.91859056e00,
        9.75344951e00,
        1.16942807e01,
        1.37138208e01,
        1.57783180e01,
        1.78549565e01,
        1.98935679e01,
        2.18502631e01,
        2.37224669e01,
        2.55036580e01,
        2.71815622e01,
        2.87508196e01,
        3.01842338e01,
        3.14663951e01,
        3.25930869e01,
    ]
    ans_sfr3 = np.array(ans_sfr3)
    d = res_sfr3 - ans_sfr3
    msg = f"{res_sfr3} {ans_sfr3} {d}"
    assert np.allclose(res_sfr3, ans_sfr3, atol=atol), msg

    res_sfr4 = sfaconc[:, 37]
    ans_sfr4 = [
        -4.26555521e-20,
        4.56008774e-02,
        2.99922765e-01,
        1.02735889e00,
        2.44958827e00,
        4.61482913e00,
        7.38275310e00,
        1.04961561e01,
        1.36797327e01,
        1.67233556e01,
        1.95047332e01,
        2.19806115e01,
        2.41479379e01,
        2.60416617e01,
        2.77139559e01,
        2.92100779e01,
        3.05620365e01,
        3.17850510e01,
        3.28886539e01,
        3.38939172e01,
        3.48144683e01,
        3.56565837e01,
        3.64265895e01,
        3.71193129e01,
        3.77329539e01,
        3.82691489e01,
    ]
    ans_sfr4 = np.array(ans_sfr4)
    d = res_sfr4 - ans_sfr4
    msg = f"{res_sfr4} {ans_sfr4} {d}"
    assert np.allclose(res_sfr4, ans_sfr4, atol=atol), msg

    # used to make results for the gwtgwt version of this problem
    # fname = os.path.join(ws, f"result_conc_lak1.txt")
    # np.savetxt(fname, res_lak1)
    # fname = os.path.join(ws, f"result_conc_sfr3.txt")
    # np.savetxt(fname, res_sfr3)
    # fname = os.path.join(ws, f"result_conc_sfr4.txt")
    # np.savetxt(fname, res_sfr4)


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
    )
    test.run()

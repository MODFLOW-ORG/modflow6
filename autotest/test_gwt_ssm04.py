"""
Test the SSM FILEINPUT option for specifying source and sink
concentrations.

Four different recharge packages are tested with the SSM FILEINPUT
1.  list-based recharge no time series
2.  array-based recharge, no time array series
3.  list-based recharge with time series
4.  array-based recharge with time array series
"""

import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["ssm04"]

nlay, nrow, ncol = 3, 5, 5
idomain_lay0 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1],
    [1, 1, 1, 1, 1],
]
idomain = np.ones((nlay, nrow, ncol), dtype=int)
idomain[0, :, :] = np.array(idomain_lay0)


def build_models(idx, test):
    perlen = [5.0, 5.0, 5.0]
    nstp = [5, 5, 5]
    tsmult = [1.0, 1.0, 1.0]
    nper = len(perlen)
    delr = 1.0
    delc = 1.0
    top = 4.0
    botm = [3.0, 2.0, 1.0]
    strt = 4.0
    hk = 1.0
    laytyp = 0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
    )

    # create iterative model solution and register the gwf model with it
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
    sim.register_ims_package(imsgwf, [gwf.name])

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

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=laytyp,
        k=hk,
        save_specific_discharge=True,
    )

    # chd files
    spd = [[(nlay - 1, nrow - 1, ncol - 1), 4.0]]
    chd = flopy.mf6.modflow.ModflowGwfchd(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="CHD-1",
    )

    # list based recharge, recharge is equal to one-based user node number
    idxrow, idxcol = np.where(idomain[0] == 1)
    recharge_rate = np.arange(nrow * ncol).reshape((nrow, ncol)) + 1
    spd = {}
    for kper in range(nper):
        rlist = []
        for i, j in zip(idxrow, idxcol):
            rlist.append([(0, i, j), recharge_rate[i, j]])
        spd[kper] = rlist
    rch1 = flopy.mf6.modflow.ModflowGwfrch(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="RCH-1",
        filename=f"{gwfname}.rch1",
    )

    # array-based rch files
    rspd = {}
    for kper in range(nper):
        rspd[kper] = recharge_rate
    rch2 = flopy.mf6.ModflowGwfrcha(
        gwf,
        print_flows=True,
        recharge=rspd,
        pname="RCH-2",
        filename=f"{gwfname}.rch2",
    )

    # list-based recharge with time series
    idxrow, idxcol = np.where(idomain[0] == 1)
    spd = []
    for i, j in zip(idxrow, idxcol):
        nodeu = i * ncol + j
        tsname = f"rch-{nodeu + 1}"
        spd.append([(0, i, j), tsname])

    tsnames = []
    for i in range(nrow):
        for j in range(ncol):
            nodeu = i * nrow + j
            tsnames.append(f"rch-{nodeu + 1}")
    ts_data = []
    totim = 0.0
    for kper in range(nper):
        totim += perlen[kper]
    for t in [0, totim]:
        ts = tuple([float(t)] + list(range(1, nrow * ncol + 1)))
        ts_data.append(ts)
    ts_dict = {
        "timeseries": ts_data,
        "time_series_namerecord": tsnames,
        "interpolation_methodrecord": [nrow * ncol * ("linear",)],
        "filename": f"{gwfname}.rch3.ts",
    }

    rch3 = flopy.mf6.modflow.ModflowGwfrch(
        gwf,
        print_flows=True,
        maxbound=len(spd),
        stress_period_data=spd,
        pname="RCH-3",
        filename=f"{gwfname}.rch3",
        timeseries=ts_dict,
    )

    # array-based rch files
    rch4 = flopy.mf6.ModflowGwfrcha(
        gwf,
        print_flows=True,
        recharge="TIMEARRAYSERIES rcharray",
        pname="RCH-4",
        filename=f"{gwfname}.rch4",
    )
    filename = f"{gwfname}.rch4.tas"
    # for now write the recharge concentration to a dat file because there
    # is a bug in flopy that will not correctly write this array as internal
    tas_array = {
        0.0: f"{gwfname}.rch4.tas.dat",
        5.0: f"{gwfname}.rch4.tas.dat",
        10.0: f"{gwfname}.rch4.tas.dat",
        15.0: f"{gwfname}.rch4.tas.dat",
    }
    time_series_namerecord = "rcharray"
    interpolation_methodrecord = "linear"
    rch4.tas.initialize(
        filename=filename,
        tas_array=tas_array,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
    )
    np.savetxt(os.path.join(ws, f"{gwfname}.rch4.tas.dat"), recharge_rate, fmt="%7.1f")

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
    )

    # create gwt model
    gwtname = "gwt_" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )
    gwt.name_file.save_flows = True

    # create iterative model solution and register the gwt model with it
    imsgwt = flopy.mf6.ModflowIms(
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

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0)

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt)

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # ssm package
    sourcerecarray = [()]
    fileinput = [
        ("RCH-1", f"{gwtname}.rch1.spc"),
        ("RCH-2", f"{gwtname}.rch2.spc"),
        ("RCH-3", f"{gwtname}.rch3.spc"),
        ("RCH-4", f"{gwtname}.rch4.spc"),
    ]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, print_flows=True, sources=sourcerecarray, fileinput=fileinput
    )

    # spc package for RCH-1
    idxrow, idxcol = np.where(idomain[0] == 1)
    recharge_concentration = np.arange(nrow * ncol).reshape((nrow, ncol)) + 1
    pd = []
    for ipos, (i, j) in enumerate(zip(idxrow, idxcol)):
        pd.append([ipos, "CONCENTRATION", recharge_concentration[i, j]])
    spc1 = flopy.mf6.ModflowUtlspc(
        gwt,
        perioddata=pd,
        maxbound=len(pd),
        filename=f"{gwtname}.rch1.spc",
    )

    # spc package for RCH-2
    idxrow, idxcol = np.where(idomain[0] == 1)
    recharge_concentration = np.arange(nrow * ncol).reshape((nrow, ncol)) + 1
    crchspd = {}
    for kper in range(nper):
        crchspd[kper] = recharge_concentration
    spc2 = flopy.mf6.ModflowUtlspca(
        gwt,
        concentration=crchspd,
        filename=f"{gwtname}.rch2.spc",
    )

    # spc package for RCH-3
    idxrow, idxcol = np.where(idomain[0] == 1)
    pd = []
    for ipos, (i, j) in enumerate(zip(idxrow, idxcol)):
        nodeu = i * ncol + j
        tsname = f"crch-{nodeu + 1}"
        pd.append([ipos, "CONCENTRATION", tsname])

    tsnames = []
    for i in range(nrow):
        for j in range(ncol):
            nodeu = i * ncol + j
            tsnames.append(f"crch-{nodeu + 1}")
    ts_data = [tuple([0.0] + list(range(1, nrow * ncol + 1)))]
    for t in [5.0, 10.0, 15.0]:
        ts = tuple([float(t)] + list(range(1, nrow * ncol + 1)))
        ts_data.append(ts)
    ts_dict = {
        "timeseries": ts_data,
        "time_series_namerecord": tsnames,
        "interpolation_methodrecord": [nrow * ncol * ("linear",)],
        "sfacrecord": [nrow * ncol * (1.0,)],
        "filename": f"{gwtname}.rch3.spc.ts",
    }
    spc3 = flopy.mf6.ModflowUtlspc(
        gwt,
        perioddata=pd,
        maxbound=len(pd),
        filename=f"{gwtname}.rch3.spc",
        timeseries=ts_dict,
        print_input=True,
    )

    # spc package for RCH-4
    spc4 = flopy.mf6.ModflowUtlspca(
        gwt,
        concentration="TIMEARRAYSERIES carray",
        filename=f"{gwtname}.rch4.spc",
        print_input=True,
    )
    filename = f"{gwtname}.rch4.spc.tas"
    # for now write the recharge concentration to a dat file because there
    # is a bug in flopy that will not correctly write this array as internal
    tas_array = {
        0.0: f"{gwtname}.rch4.spc.tas.dat",
        5.0: f"{gwtname}.rch4.spc.tas.dat",
        10.0: f"{gwtname}.rch4.spc.tas.dat",
        15.0: f"{gwtname}.rch4.spc.tas.dat",
    }
    time_series_namerecord = "carray"
    interpolation_methodrecord = "linear"
    spc4.tas.initialize(
        filename=filename,
        tas_array=tas_array,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
    )
    recharge_concentration = np.arange(nrow * ncol).reshape((nrow, ncol)) + 1
    np.savetxt(
        os.path.join(ws, f"{gwtname}.rch4.spc.tas.dat"),
        recharge_concentration,
        fmt="%7.1f",
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )

    obs_data = {
        f"{gwtname}.obs.csv": [
            ("(1-1-1)", "CONCENTRATION", (0, 0, 0)),
            ("(1-5-5)", "CONCENTRATION", (nlay - 1, nrow - 1, ncol - 1)),
        ],
    }

    obs_package = flopy.mf6.ModflowUtlobs(
        gwt,
        pname=f"{gwtname}.obs",
        digits=10,
        print_input=True,
        continuous=obs_data,
    )

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    name = test.name
    gwtname = "gwt_" + name

    # load concentration file
    fpth = os.path.join(test.workspace, f"{gwtname}.ucn")
    cobj = flopy.utils.HeadFile(fpth, precision="double", text="CONCENTRATION")
    conc = cobj.get_data()

    # load transport budget file
    fpth = os.path.join(test.workspace, f"{gwtname}.cbc")
    bobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    ssmbudall = bobj.get_data(text="SOURCE-SINK MIX")
    times = cobj.get_times()

    print(times)
    for itime, totim in enumerate(times):
        print(f"Checking records for time {totim}")

        # Check records for each of the four recharge packages
        ssmbud = ssmbudall[itime]
        istart = 0
        for irchpak in [1, 2, 3, 4]:
            print(f"  Checking records for recharge package {irchpak}")
            istop = istart + 23

            print(ssmbud[istart:istop])

            print("    Checking id1")
            id1 = ssmbud[istart:istop]["node"]
            id1a = [
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                14,
                15,
                16,
                17,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
            ]
            assert np.allclose(id1, id1a), f"{id1} /= {id1a}"

            print("    Checking id2")
            id2 = ssmbud[istart:istop]["node2"]
            if irchpak in [1, 3]:
                # recharge packages 1 and 3 are list-based with 23 entries
                id2a = np.arange(23) + 1
            elif irchpak in [2, 4]:
                # recharge packages 2 and 4 are array-based with 25 entries
                id2a = id1a
            assert np.allclose(id2, id2a), f"q: {id2} /= {id2a}"

            print(f"    Checking q for irchpak {irchpak}")
            q = ssmbud[istart:istop]["q"]
            if irchpak in [2, 3]:
                qa = [float(a) ** 2 for a in id1]
            else:
                qa = [float(a) ** 2 for a in id1a]
            for i in range(23):
                print(f"{i + 1} {id1[i]} {id2[i]} {q[i]} {qa[i]}")
            assert np.allclose(q, qa), f"q: {q} /=\n {qa}"

            istart = istop


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

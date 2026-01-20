"""
Test that zonebudget works on a cell budget file from GWT
https://github.com/MODFLOW-ORG/modflow6/discussions/1181
"""

import os
from pathlib import Path

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["zbud6_zb01"]
htol = None
dtol = 1e-3
budtol = 1e-2
bud_lst = ["STORAGE-AQUEOUS_IN", "STORAGE-AQUEOUS_OUT"]
zone_lst = []
for n in bud_lst:
    s = n.replace("_", "-")
    zone_lst.append(s)


nlay, nrow, ncol = 5, 10, 20
nper = 1
delr = 1.0
delc = 1.0
delz = 1.0
top = 1.0
botm = np.linspace(top - delz, top - nlay * delz, nlay)
strt = 1.0
hk = 1.0
laytyp = 0
porosity = 0.1
qwell = 1.0
specific_discharge = qwell / delr / delz
timetoend = float(ncol) * delc * porosity / specific_discharge
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol


def build_models(idx, test):
    name = cases[idx]
    perlen = [timetoend]
    nstp = [50]
    tsmult = [1.0]
    steady = [True]

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-6, 1.0

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
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        model_nam_file=f"{gwfname}.nam",
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
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    # chd and wel info
    chdlist = []
    wellist = []
    for k in range(nlay):
        for i in range(nrow):
            for j in range(ncol):
                if j == ncol - 1:
                    chdlist.append([(k, i, j), 0.0])
                if j == 0:
                    wellist.append([(k, i, j), qwell, 1.0])
    c = {0: chdlist}
    w = {0: wellist}

    # grid discretization
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwfname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{gwfname}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        icelltype=laytyp,
        xt3doptions=[()],
        k=hk,
        k33=hk,
        save_specific_discharge=True,
    )

    # chd package
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(c),
        stress_period_data=c,
        save_flows=False,
        pname="CHD-1",
    )

    # wel package
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=len(w),
        stress_period_data=w,
        save_flows=False,
        auxiliary="CONCENTRATION",
        pname="WEL-1",
    )

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

    # gwt grid discretization
    dis = flopy.mf6.ModflowGwtdis(
        gwt,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        idomain=np.ones((nlay, nrow, ncol), dtype=int),
        filename=f"{gwtname}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwtic(gwt, strt=0.0, filename=f"{gwtname}.ic")

    # advection
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme="upstream", filename=f"{gwtname}.adv")

    # mass storage and transfer
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=0.1)

    # sources
    sourcerecarray = [("WEL-1", "AUX", "CONCENTRATION")]
    ssm = flopy.mf6.ModflowGwtssm(
        gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm"
    )

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "LAST")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
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
    ws = Path(test.workspace)

    # build zonebudget files
    # start with 1 since budget isn't calculated for zone 0
    zones = [k + 1 for k in range(nlay)]
    nzones = len(zones)
    with open(ws / "zonebudget.nam", "w") as f:
        f.write("BEGIN ZONEBUDGET\n")
        f.write(f"  BUD gwt_{test.name}.cbc\n")
        f.write(f"  ZON {test.name}.zon\n")
        f.write(f"  GRB gwf_{test.name}.dis.grb\n")
        f.write("END ZONEBUDGET\n")

    with open(ws / f"{test.name}.zon", "w") as f:
        f.write("BEGIN DIMENSIONS\n")
        f.write(f"  NCELLS {size3d}\n")
        f.write("END DIMENSIONS\n\n")
        f.write("BEGIN GRIDDATA\n")
        f.write("  IZONE LAYERED\n")
        for k in range(nlay):
            f.write(f"    CONSTANT {zones[k]:>10d}\n")
        f.write("END GRIDDATA\n")

    # run zonebudget
    success, buff = flopy.run_model(
        test.targets["zbud6"],
        "zonebudget.nam",
        model_ws=ws,
        silent=False,
        report=True,
    )

    assert success
    test.success = success

    # read data from csv file
    zbd = np.genfromtxt(
        ws / "zonebudget.csv", names=True, delimiter=",", deletechars=""
    )

    # sum the data for all zones
    nentries = int(zbd.shape[0] / nzones)
    zbsum = np.zeros(nentries, dtype=zbd.dtype)
    static = ["totim", "kstp", "kper"]
    ipos = 0
    ion = 0
    for t in zbd:
        for name in zbd.dtype.names:
            if name in static:
                zbsum[name][ipos] = t[name]
            elif name == "zone":
                zbsum[name][ipos] = 0
            else:
                zbsum[name][ipos] += t[name]
        ion += 1
        if ion == nzones:
            ipos += 1
            ion = 0

    # get results from listing file
    # todo: should flopy have a subclass for GWT list file?
    budl = flopy.utils.mflistfile.ListBudget(
        ws / f"gwt_{os.path.basename(test.name)}.lst",
        budgetkey="MASS BUDGET FOR ENTIRE MODEL",
    )
    names = list(bud_lst)
    found_names = budl.get_record_names()
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]

    # get results from cbc file
    cbc_bud = ["STORAGE-AQUEOUS"]
    d = np.recarray(nbud, dtype=dtype)
    for key in bud_lst:
        d[key] = 0.0
    cobj = flopy.utils.CellBudgetFile(
        ws / f"gwt_{os.path.basename(test.name)}.cbc", precision="double"
    )
    rec = cobj.list_records()
    kk = cobj.get_kstpkper()
    times = cobj.get_times()
    for i, (k, t) in enumerate(zip(kk, times)):
        for text in cbc_bud:
            qin = 0.0
            qout = 0.0
            v = cobj.get_data(kstpkper=k, text=text)[0]
            if isinstance(v, np.recarray):
                vt = np.zeros(size3d, dtype=float)
                for jdx, node in enumerate(v["node"]):
                    vt[node - 1] += v["q"][jdx]
                v = vt.reshape(shape3d)
            for kk in range(v.shape[0]):
                for ii in range(v.shape[1]):
                    for jj in range(v.shape[2]):
                        vv = v[kk, ii, jj]
                        if vv < 0.0:
                            qout -= vv
                        else:
                            qin += vv
            d["totim"][i] = t
            d["time_step"][i] = k[0]
            d["stress_period"] = k[1]
            key = f"{text}_IN"
            d[key][i] = qin
            key = f"{text}_OUT"
            d[key][i] = qout

    # calculate absolute difference
    diff = np.zeros((nbud, len(bud_lst)), dtype=float)
    for i, key in enumerate(bud_lst):
        diff[:, i] = d0[key] - d[key]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-budget difference ({diffmax}) "

    # write summary
    with open(ws / f"{os.path.basename(test.name)}.bud.cmp.out", "w") as f:
        for i in range(diff.shape[0]):
            if i == 0:
                line = f"{'TIME':>10s}"
                for key in bud_lst:
                    line += f"{key + '_LST':>25s}"
                    line += f"{key + '_CBC':>25s}"
                    line += f"{key + '_DIF':>25s}"
                f.write(line + "\n")
            line = f"{d['totim'][i]:10g}"
            for ii, key in enumerate(bud_lst):
                line += f"{d0[key][i]:25g}"
                line += f"{d[key][i]:25g}"
                line += f"{diff[i, ii]:25g}"
            f.write(line + "\n")

    # compare zone budget output to cbc output
    diffzb = np.zeros((nbud, len(bud_lst)), dtype=float)
    for i, (key0, key) in enumerate(zip(zone_lst, bud_lst)):
        diffzb[:, i] = zbsum[key0] - d[key]
    diffzbmax = np.abs(diffzb).max()
    msg += f"\nmaximum absolute zonebudget-cell by cell difference ({diffzbmax}) "

    # write summary
    with open(ws / f"{os.path.basename(test.name)}.zbud.cmp.out", "w") as f:
        for i in range(diff.shape[0]):
            if i == 0:
                line = f"{'TIME':>10s}"
                for key in bud_lst:
                    line += f"{key + '_ZBUD':>25s}"
                    line += f"{key + '_CBC':>25s}"
                    line += f"{key + '_DIF':>25s}"
                f.write(line + "\n")
            line = f"{d['totim'][i]:10g}"
            for ii, (key0, key) in enumerate(zip(zone_lst, bud_lst)):
                line += f"{zbsum[key0][i]:25g}"
                line += f"{d[key][i]:25g}"
                line += f"{diffzb[i, ii]:25g}"
            f.write(line + "\n")

    if diffmax > budtol or diffzbmax > budtol:
        test.success = False
        msg += f"\n...exceeds {budtol}"
        assert diffmax < budtol and diffzbmax < budtol, msg
    else:
        test.success = True
        print("    " + msg)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        htol=htol,
    )
    test.run()

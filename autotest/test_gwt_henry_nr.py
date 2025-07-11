"""
This is the Henry, Newton-Raphson problem described by Langevin et al (2020)
with a 20 by 40 grid instead of the 40 by 80 grid described in the paper.
There is freshwater inflow on the left and a sloping sea boundary on the
right with moves up and down according to a simple sine function.  GHBs
and DRNs alternate and move up and down along the boundary to represent
the effects of tides on the aquifer.
"""

import os

import flopy
import numpy as np
import pytest
from framework import DNODATA, TestFramework

cases = ["henrynr01"]

# global model variables
nlay = 20
ncol = 40
fx = 0.5
fz = 0.5
lx = 2
lz = 1
sealevel = 0.85
amplitude = 0.14
frequency = 4
wellfact = 0.25


def get_idomain(nlay, nrow, ncol, lx, lz, fx, fz):
    idomain = np.ones((nlay, nrow, ncol), dtype=int)
    x1 = fx * lx
    y1 = lz
    x2 = lx
    y2 = fz * lz
    slope = (y2 - y1) / (x2 - x1)
    b = y1 - slope * x1

    delr = lx / ncol
    delv = lz / nlay
    xcenters = np.linspace(delr / 2, lx - delr / 2, ncol)
    zcenters = np.linspace(lz - delv / 2, delv / 2, nlay)

    for k in range(nlay):
        zc = zcenters[k]
        for j in range(ncol):
            xc = xcenters[j]
            zedge = slope * xc + b
            if zc > zedge:
                idomain[k, 0, j] = 0

    kidm0, iidmn0, jidmn0 = np.where(idomain == 0)
    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k, 0, j - 1] == 1:
            idomain[k, 0, j - 1] = 2

    for k, j in zip(kidm0, jidmn0):
        if idomain[k, 0, j] == 0 and idomain[k + 1, 0, j] == 1:
            idomain[k + 1, 0, j] = 3

    return idomain


def sinfunc(a, b, c, d, x):
    return a * np.sin(b * (x - c)) + d


def get_model(ws, name, array_input=False):
    nrow = 1
    delr = lx / ncol
    delc = 1.0
    top = lz
    delz = lz / nlay
    botm = list(top - np.arange(delz, nlay * delz + delz, delz))

    perlen = [0.25] + 1000 * [0.001]
    nper = len(perlen)
    nstp = [250] + 1000 * [1]
    tsmult = 1.0
    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp[i], tsmult))

    nouter, ninner = 200, 50
    hclose, rclose, relax = 1e-7, 1e-5, 0.97

    # build MODFLOW 6 files
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    sim.name_file.continue_ = False

    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name

    # --------------------  FLOW --------------------

    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, newtonoptions="NEWTON")

    imsgwf = flopy.mf6.ModflowIms(
        sim,
        print_option="summary",
        csv_outer_output_filerecord=gwfname + ".ims.csv",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename=f"{gwfname}.ims",
        no_ptcrecord=True,
    )
    sim.register_ims_package(imsgwf, [gwfname])

    idomain = get_idomain(nlay, nrow, ncol, lx, lz, fx=fx, fz=fz)
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
    ic = flopy.mf6.ModflowGwfic(gwf, strt=lz)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        xt3doptions=False,
        save_flows=True,
        save_specific_discharge=True,
        icelltype=1,
        k=864.0,
    )

    sto = flopy.mf6.ModflowGwfsto(
        gwf, sy=0.35, iconvert=1, steady_state=[False], transient=[True]
    )

    pd = [(0, 0.7, 0.0, gwtname, "none")]
    buy = flopy.mf6.ModflowGwfbuy(gwf, hhformulation_rhs=False, packagedata=pd)

    # drn and ghb
    kidx, iidx, jidx = np.where(idomain > 1)
    xcellcenters = gwf.modelgrid.xcellcenters
    zcellcenters = gwf.modelgrid.zcellcenters
    botm = dis.botm.get_data()
    dt = 0.001
    times = np.arange(dt, 1.0 + dt, dt)
    sealevelts = [sealevel] + list(
        sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, times)
    )
    drnspd = {}
    if array_input:
        bheadspd = {}
        ghbcondspd = {}
        ghbauxspd = {}
        elevspd = {}
        drncondspd = {}
        drnauxspd = {}
    else:
        ghbspd = {}
        drnspd = {}
    for kper in range(nper):
        if kper == 0:
            sl = sealevel
        else:
            sl = sealevelts[kper]
        sl = np.round(sl, decimals=8)
        if array_input:
            bhead = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            ghbcond = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            ghbconc = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            ghbdens = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            elev = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            drncond = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
            drnconc = np.full((nlay, nrow, ncol), DNODATA, dtype=float)
        else:
            ghblist = []
            drnlist = []
        ghbbnd = 0
        drnbnd = 0
        for k, i, j in zip(kidx, iidx, jidx):
            zcell = zcellcenters[k, i, j]
            cond = 864.0 * (delz * delc) / (0.5 * delr)
            if zcell > sl:
                if array_input:
                    elev[k, i, j] = zcell
                    drncond[k, i, j] = 864.0
                    drnconc[k, i, j] = 0.0
                else:
                    drnlist.append([(k, i, j), zcell, 864.0, 0.0])
                drnbnd += 1
            else:
                if array_input:
                    bhead[k, i, j] = sl
                    ghbcond[k, i, j] = 864.0
                    ghbconc[k, i, j] = 35.0
                    ghbdens[k, i, j] = 1024.5
                else:
                    ghblist.append([(k, i, j), sl, 864.0, 35.0, 1024.5])
                ghbbnd += 1
        if ghbbnd > 0:
            if array_input:
                bheadspd[kper] = bhead
                ghbcondspd[kper] = ghbcond
                ghbauxspd[kper] = [ghbconc, ghbdens]
            else:
                ghbspd[kper] = ghblist
        if drnbnd > 0:
            if array_input:
                elevspd[kper] = elev
                drncondspd[kper] = drncond
                drnauxspd[kper] = [drnconc]
            else:
                drnspd[kper] = drnlist

    # drn
    if array_input:
        drn1 = flopy.mf6.ModflowGwfdrng(
            gwf,
            print_input=True,
            print_flows=True,
            save_flows=False,
            maxbound=11,
            pname="DRN-1",
            auxiliary=["CONCENTRATION"],
            elev=elevspd,
            cond=drncondspd,
            aux=drnauxspd,
        )
    else:
        drn1 = flopy.mf6.ModflowGwfdrn(
            gwf,
            stress_period_data=drnspd,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="DRN-1",
            auxiliary="CONCENTRATION",
        )

    # ghb / ghbg
    if array_input:
        ghb1 = flopy.mf6.ModflowGwfghbg(
            gwf,
            print_input=True,
            print_flows=True,
            save_flows=False,
            maxbound=20,
            pname="GHB-1",
            auxiliary=["CONCENTRATION", "DENSITY"],
            bhead=bheadspd,
            cond=ghbcondspd,
            aux=ghbauxspd,
        )
    else:
        ghb1 = flopy.mf6.ModflowGwfghb(
            gwf,
            stress_period_data=ghbspd,
            print_input=True,
            print_flows=True,
            save_flows=False,
            pname="GHB-1",
            auxiliary=["CONCENTRATION", "DENSITY"],
        )

    wellist1 = []
    qwell = 5.7024 * wellfact
    qwell = qwell / nlay
    for k in range(nlay):
        wellist1.append([(k, 0, 0), qwell, 0.0])
    wel1 = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=wellist1,
        print_input=True,
        print_flows=True,
        save_flows=False,
        pname="WEL-1",
        auxiliary="CONCENTRATION",
        filename=f"{gwfname}.wel",
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    # --------------------  TRANSPORT --------------------

    # create gwt model
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="all",
        outer_dvclose=hclose,
        outer_maximum=100,
        inner_maximum=50,
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
    ic = flopy.mf6.ModflowGwtic(gwt, strt=35.0)

    # advection
    scheme = "UPSTREAM"
    # scheme = 'TVD'
    adv = flopy.mf6.ModflowGwtadv(gwt, scheme=scheme)

    # dispersion
    diffusion_only = True
    if diffusion_only:
        diffc = 0.57024
        alh = 0.0
        ath = 0.0
        xt3d = False
    else:
        diffc = 0.0
        alh = 0.1
        ath = 0.01
        xt3d = True
    xt3d_off = not xt3d
    dsp = flopy.mf6.ModflowGwtdsp(
        gwt, xt3d_off=xt3d_off, diffc=diffc, alh=alh, ath1=ath
    )

    # mass storage and transfer
    porosity = 0.35
    mst = flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)

    # sources
    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("WEL-1", "AUX", "CONCENTRATION"),
    ]
    ssm = flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)

    # fmi
    # this problem does not want to converge unless the flow correction
    # is implemented.  It looks like very small flow errors can wreak
    # havoc on the transport solution.  The correction seems to mitigate
    # these problems and at least keep the solution stable
    fmi = flopy.mf6.ModflowGwtfmi(gwt, flow_imbalance_correction=True)

    # output control
    oc = flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "ALL")],
    )

    # --------------------  EXCHANGE --------------------

    # GWF GWT exchange
    gwfgwt = flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=gwfname,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim


def get_patch_collection(modelgrid, head, conc, cmap="jet", zorder=None):
    # create patches for each cell
    import matplotlib.collections
    import matplotlib.patches

    xv, yv, zv = modelgrid.xyzvertices
    botm = modelgrid.botm
    patches = []
    for k in range(modelgrid.nlay):
        for j in range(modelgrid.ncol):
            x0 = xv[0, j]
            x1 = xv[0, j + 1]
            z0 = zv[k, 0, j]
            z0 = min(z0, head[k, 0, j])
            z0 = max(z0, botm[k, 0, j])
            z1 = zv[k + 1, 0, j]
            poly = [[x0, z0], [x1, z0], [x1, z1], [x0, z1], [x0, z0]]
            patch = matplotlib.patches.Polygon(
                poly, closed=True, edgecolor="k", facecolor="red"
            )
            patches.append(patch)
    pc = matplotlib.collections.PatchCollection(patches, cmap=cmap, zorder=zorder)
    pc.set_array(conc.flatten())
    return pc


def plot_output(idx, test):
    ws = test.workspace
    name = test.name
    sim = test.sims[0]
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)
    headall = gwf.output.head().get_alldata()
    concall = gwt.output.concentration().get_alldata()

    import matplotlib.patches
    import matplotlib.pyplot as plt

    # sea level figure
    fig = plt.figure(figsize=(5, 5), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    dt = 0.001
    times = np.arange(dt, 1.0 + dt, dt)
    sealevelts = 250 * [sealevel] + list(
        sinfunc(amplitude, frequency * 2 * np.pi, 0, sealevel, times)
    )
    simtime = np.linspace(dt, 1.25 + dt, 1250)
    plt.plot(simtime, sealevelts)
    ax.set_xlabel("TIME, IN DAYS", fontsize=6)
    ax.set_ylabel("SEAWATER LEVEL, IN METERS", fontsize=6)
    fname = "fig-sealevel.pdf"
    fname = os.path.join(ws, fname)
    plt.savefig(fname, bbox_inches="tight")

    # results plot
    fig = plt.figure(figsize=(8.5, 11.0), dpi=300)
    botm = gwf.modelgrid.botm
    levels = [0.01, 0.1, 0.5, 0.9, 0.99]
    times2plot = [249, np.argmax(sealevelts), np.argmin(sealevelts)]
    nplots = len(times2plot)
    figtxt = ["(a)", "(b)", "(c)"]
    for ifig, itime in enumerate(times2plot):
        ax = fig.add_subplot(nplots, 1, ifig + 1, aspect="equal")
        conc = concall[itime]
        head = headall[itime]
        conc = np.ma.masked_greater(conc, 1e20)
        conc = np.ma.masked_where(head < botm, conc)
        sl = sealevelts[itime]
        # sea polygon
        seapoly = np.array([[lx * fx, sl], [lx, sl], [lx, 0]])
        patch = matplotlib.patches.Polygon(
            seapoly, closed=True, facecolor="darkred", zorder=0
        )
        ax.add_patch(patch)
        # aquifer polygon
        aqpoly = np.array([[0, 0], [lx, 0], [lx, fz * lz], [lx * fx, lz], [0, lz]])
        patch = matplotlib.patches.Polygon(
            aqpoly, closed=True, facecolor=".7", zorder=1
        )
        ax.add_patch(patch)
        # model cross section
        xs = flopy.plot.PlotCrossSection(gwf, line={"row": 0}, ax=ax)
        # color filled model cells
        pc = get_patch_collection(gwf.modelgrid, head, conc, zorder=2)
        pc.set_clim(0, 35.0)
        ax.add_collection(pc)
        # model grid
        xs.plot_grid(linewidths=0.5)
        # concentration contours
        cs = plt.contour(
            np.flipud(conc[:, 0, :] / 35.0),
            extent=[0, lx, 0, lz],
            levels=levels,
            colors="white",
            zorder=2,
        )
        ax.clabel(cs, fontsize=6, fmt="%1.2f")
        # labels and title
        if ifig == nplots - 1:
            ax.set_xlabel("DISTANCE, IN METERS", fontsize=6)
        ax.set_ylabel("ELEVATION, IN METERS", fontsize=6)
        ttl = f"TIME = {simtime[itime]:.3f} days"
        ax.set_title(ttl, fontsize=6)
        ax.text(1.9, 1.025, figtxt[ifig], fontsize=6)

    fname = "fig-concplots.pdf"
    fname = os.path.join(ws, fname)
    plt.savefig(fname, bbox_inches="tight")


def check_output(ws, name, sim):
    gwfname = "gwf_" + name
    gwtname = "gwt_" + name
    gwf = sim.get_model(gwfname)
    gwt = sim.get_model(gwtname)
    headobj = gwf.output.head()
    concobj = gwt.output.concentration()

    # extract 10 simulated heads and concs for cell (0, 0, 20)
    hsim = headobj.get_ts((0, 0, 20))[::125, 1]
    csim = concobj.get_ts((0, 0, 20))[::125, 1]
    hans = np.array(
        [
            0.88909027,
            0.85382205,
            0.85710524,
            0.85227865,
            0.85692375,
            0.85245585,
            0.85698074,
            0.8524655,
            0.85698465,
            0.85246506,
        ]
    )
    cans = np.array(
        [
            35.0,
            15.25952065,
            8.98272748,
            9.12073523,
            7.73799374,
            7.56138739,
            7.37415009,
            7.46278052,
            7.34050542,
            7.49958207,
        ]
    )

    errmsg = f"heads not right for cell (0, 0, 20):\n{hsim}\n{hans}"
    assert np.allclose(hsim, hans, atol=1.0e-3), errmsg

    errmsg = f"concs not right for cell (0, 0, 20):\n{csim}\n{cans}"
    assert np.allclose(hsim, hans, atol=1.0e-3), errmsg


def build_models(idx, test):
    # build MODFLOW 6 files
    ws = test.workspace
    name = cases[idx]
    sim = get_model(ws, name)

    # build comparison array_input model
    ws = os.path.join(test.workspace, "mf6")
    mc = get_model(ws, name, array_input=True)

    return sim, mc


def check_outputs(idx, test):
    name = cases[idx]
    sim = test.sims[0]

    # check output MODFLOW 6 files
    ws = test.workspace
    check_output(ws, name, sim)

    # check output comparison array_input model
    ws = os.path.join(test.workspace, "mf6")
    check_output(ws, name, sim)


@pytest.mark.slow
@pytest.mark.developmode
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_outputs(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        compare="mf6",
    )
    test.run()

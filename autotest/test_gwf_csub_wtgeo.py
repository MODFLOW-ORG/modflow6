import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = [
    "csub_wtgeoa",
    "csub_wtgeob",
    "csub_wtgeoc",
    "csub_wtgeod",
    "csub_wtgeoe",
    "csub_wtgeof",
    "csub_wtgeog",
]
constantcv = [True for _ in range(len(cases))]
cmppth = "mf6_regression"
compare = [True for _ in range(len(cases))]
tops = [0.0, 0.0, 150.0, 0.0, 0.0, 150.0, 150.0]
ump = [None, None, True, None, True, None, True]
iump = [0, 0, 1, 0, 1, 0, 1]
eslag = [True for _ in range(len(cases) - 2)] + 2 * [False]
# eslag = [True, True, True, False, True, False, False]
headformulation = [True, False, False, True, True, False, False]
ndc = [None, None, None, 19, 19, 19, 19]
delay = [False, False, False, True, True, True, True]
# newton = ["", "", "", "", "", None, ""]
newton = ["NEWTON" for _ in range(len(cases))]

htol = [None, None, None, 0.2, None, None, None]
dtol = 1e-3
budtol = 1e-2
paktest = "csub"

# static model data
# temporal discretization
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# spatial discretization data
nlay, nrow, ncol = 3, 10, 10
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1000.0, 2000.0
botm = [-100, -150.0, -350.0]
strt = 0.0
hnoflo = 1e30
hdry = -1e30

# calculate hk
hk1fact = 1.0 / 50.0
hk1 = np.ones((nrow, ncol), dtype=float) * 0.5 * hk1fact
hk1[0, :] = 1000.0 * hk1fact
hk1[-1, :] = 1000.0 * hk1fact
hk1[:, 0] = 1000.0 * hk1fact
hk1[:, -1] = 1000.0 * hk1fact
hk = [20.0, hk1, 5.0]

# calculate vka
vka = [1e6, 7.5e-5, 1e6]

# all cells are active and layer 1 is convertible
ib = 1

# solver options
nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0
# newtonoptions = None
imsla = "BICGSTAB"

# chd data
c = []
c6 = []
ccol = [3, 4, 5, 6]
for j in ccol:
    c.append([0, nrow - 1, j, strt, strt])
    c6.append([(0, nrow - 1, j), strt])
cd = {0: c}
cd6 = {0: c6}
maxchd = len(cd[0])

# pumping well data
wr = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3]
wc = [0, 1, 8, 9, 0, 9, 0, 9, 0, 0]
wrp = [2, 2, 3, 3]
wcp = [5, 6, 5, 6]
wq = [-14000.0, -8000.0, -5000.0, -3000.0]
d = []
d6 = []
for r, c, q in zip(wrp, wcp, wq):
    d.append([2, r, c, q])
    d6.append([(2, r, c), q])
wd = {1: d}
wd6 = {1: d6}
maxwel = len(wd[1])

# recharge data
q = 3000.0 / (delr * delc)
v = np.zeros((nrow, ncol), dtype=float)
for r, c in zip(wr, wc):
    v[r, c] = q
rech = {0: v}

# storage and compaction data
sgm = 1.7
sgs = 2.0
void = 0.82
preconhead = -7.0
theta = void / (1.0 + void)
sw = 4.65120000e-10 * 9806.65000000 * theta
sy = [0.1, 0.0, 0.0]
ske = [6e-6, 3e-6, 6e-6]
skv = [6e-4, 3e-4, 6e-4]
cg_ske_cr = [ske[0], 0, ske[-1]]
kv = 1e-6

bdb = [45.0, 0, 90.0]
facib = [0.6, 1.0, 0.6]
facsk = [0.4, 0.0, 0.4]
# facndb = [0.15, 1., 0.15]
# facdb = [0.45, 0., 0.45]

dp = [[kv, ske[0], skv[0]]]
rnb = [7.635, 0.0, 17.718]
dhc = [preconhead for n in range(nlay)]
dstart = [strt for n in range(nlay)]
dz = [5.894, 0.0, 5.08]
nz = [1, 0, 1]


# calculate geostatic and effective stress
def calc_stress(sgm0, sgs0, h, bt):
    geo = []
    for k in range(nlay):
        top = bt[k]
        bot = bt[k + 1]
        ht = h
        if ht > top:
            gs = (top - bot) * sgs0
        elif ht < bot:
            gs = (top - bot) * sgm0
        else:
            gs = ((top - ht) * sgm0) + ((ht - bot) * sgs0)
        geo.append(gs)
    # calculate total geostatic stress at bottom of layer
    for k in range(1, nlay):
        geo[k] += geo[k - 1]
    # calculate effective stress at the bottom of the layer
    es = []
    for k in range(nlay):
        es.append(geo[k] - (h - bt[k + 1]))
    return geo, es


def get_model(idx, ws):
    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration=imsla,
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
    )

    # create gwf model
    top = tops[idx]
    zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
    elevs = [top] + botm

    if top == 0:
        laytyp = [0, 0, 0]
    else:
        laytyp = [1, 0, 0]

    # calculate sk, ndb, and db factors
    facdb = []
    facndb = []
    for k in range(nlay):
        bt = zthick[k]
        f = bdb[k] / bt
        facdb.append(f)
        bnd = bt * (facib[k] - f)
        f = bnd / bt
        facndb.append(f)

    # csub packagedata container counter
    if headformulation[idx]:
        head_based = True
        sgmt = None
        sgst = None
    else:
        head_based = None
        sgmt = sgm
        sgst = sgs

    # fill preconsolidation stress with preconsolidation head
    # calculate preconsolidation stress, if necessary
    pcs = [preconhead for k in range(nlay)]
    gs, es = calc_stress(sgm, sgs, preconhead, elevs)
    if not headformulation[idx]:
        pcs = es

    # create no delay bed packagedata entries
    sub6 = []
    ibcno = 0
    nndb = 0
    ln = []
    sfv = []
    sfe = []
    hc = []
    thickib0 = []
    cdelays = "nodelay"
    for k in range(nlay):
        b = zthick[k] * facndb[k]
        if b <= 0.0:
            continue
        nndb += 1
        ln.append(k)
        thickib0.append(b)
        if headformulation[idx]:
            sfv.append(skv[k] * b)
            sfe.append(ske[k] * b)
            hc.append(pcs[k])
        else:
            sfv.append(skv[k])
            sfe.append(ske[k])
        for i in range(nrow):
            for j in range(ncol):
                # skip constant head cells
                if k == 0 and i == nrow - 1 and j in ccol:
                    continue
                # create nodelay entry
                # no delay beds
                d = [
                    ibcno,
                    (k, i, j),
                    cdelays,
                    pcs[k],
                    b,
                    1.0,
                    skv[k],
                    ske[k],
                    theta,
                    999.0,
                    -999.0,
                ]
                sub6.append(d)
                ibcno += 1

    if delay[idx]:
        cdelays = "delay"
    else:
        cdelays = "nodelay"
    ndb = 0
    nmz = 0
    ldn = []
    rnbsub = []
    nzsub = []
    dzsub = []
    dhcsub = []
    dstartsub = []
    for k in range(nlay):
        b = zthick[k] * facdb[k]
        if b <= 0.0:
            continue
        if delay[idx]:
            ndb += 1
            nmz = 1
            bb = dz[k]
            rnbt = rnb[k]
            rnbsub.append(rnbt)
            dzsub.append(dz[k])
            dhcsub.append(dhc[k])
            dstartsub.append(dstart[k])
            hib = strt
            nzsub.append(1)
            ldn.append(k)
        else:
            bb = b
            rnbt = 1
            hib = -999.0
            nndb += 1
            ln.append(k)
            thickib0.append(b)
            if headformulation[idx]:
                sfv.append(skv[k] * b)
                sfe.append(ske[k] * b)
                hc.append(pcs[k])
            else:
                sfv.append(skv[k])
                sfe.append(ske[k])
        for i in range(nrow):
            for j in range(ncol):
                # skip constant head cells
                if k == 0 and i == nrow - 1 and j in ccol:
                    continue
                # create nodelay entry
                # no delay beds
                d = [
                    ibcno,
                    (k, i, j),
                    cdelays,
                    pcs[k],
                    bb,
                    rnbt,
                    skv[k],
                    ske[k],
                    theta,
                    kv,
                    hib,
                ]
                sub6.append(d)
                ibcno += 1

    # add coarse-grained component
    for k in range(nlay):
        b = zthick[k] * facsk[k]
        if b <= 0.0:
            continue
        nndb += 1
        ln.append(k)
        thickib0.append(b)
        if headformulation[idx]:
            sfv.append(ske[k] * b)
            sfe.append(ske[k] * b)
            hc.append(pcs[k])
        else:
            sfv.append(ske[k])
            sfe.append(ske[k])

    maxcsub = len(sub6)

    # water compressibility cannot be compared for cases where the material
    # properties are adjusted since the porosity changes in mf6
    if iump[idx] == 0:
        beta = 4.6512e-10
        wc = sw
    else:
        beta = 0.0
        wc = 0.0

    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions=newton[idx])

    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=False,
        # dev_modflowusg_upstream_weighted_saturation=True,
        icelltype=laytyp,
        k=hk,
        k33=vka,
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        storagecoefficient=True,
        steady_state={0: True},
        transient={1: True},
    )

    # recharge
    rch = flopy.mf6.ModflowGwfrcha(gwf, readasarrays=True, recharge=rech)

    # wel file
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=maxwel,
        stress_period_data=wd6,
        save_flows=False,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
    )
    # ibc files
    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        specified_initial_interbed_state=True,
        ndelaycells=ndc[idx],
        head_based=head_based,
        update_material_properties=ump[idx],
        effective_stress_lag=eslag[idx],
        save_flows=True,
        ninterbeds=maxcsub,
        sgm=sgmt,
        sgs=sgst,
        cg_theta=theta,
        cg_ske_cr=cg_ske_cr,
        beta=beta,
        packagedata=sub6,
    )
    obspos = [(0, 4, 4), (1, 4, 4), (2, 4, 4)]
    obstype = [
        "compaction-cell",
        "gstress-cell",
        "estress-cell",
        "ske-cell",
        "sk-cell",
        "csub-cell",
    ]
    obstag = ["tcomp", "gs", "es", "ske", "sk", "csub"]
    obsarr = []
    for iobs, cobs in enumerate(obstype):
        for jobs, otup in enumerate(obspos):
            otag = f"{obstag[iobs]}{jobs + 1}"
            obsarr.append((otag, cobs, otup))

    orecarray = {}
    orecarray["csub_obs.csv"] = obsarr
    csub_obs_package = csub.obs.initialize(
        filename=opth, digits=10, print_input=True, continuous=orecarray
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)  # modflow6 files
    mc = get_model(idx, os.path.join(test.workspace, cmppth))  # build comparison files
    return sim, mc


def check_output(idx, test):
    if compare[idx]:
        # MODFLOW 6 total compaction results
        fpth = os.path.join(test.workspace, "csub_obs.csv")
        try:
            tc = np.genfromtxt(fpth, names=True, delimiter=",")
        except:
            assert False, f'could not load data from "{fpth}"'

        # comparison total compaction results
        fpth = os.path.join(test.workspace, cmppth, "csub_obs.csv")
        try:
            tc0 = np.genfromtxt(fpth, names=True, delimiter=",")
        except:
            assert False, f'could not load data from "{fpth}"'

        # calculate maximum absolute error
        diff = tc["TCOMP3"] - tc0["TCOMP3"]
        diffmax = np.abs(diff).max()
        msg = f"maximum absolute total-compaction difference ({diffmax}) "

        # write summary
        fpth = os.path.join(
            test.workspace, f"{os.path.basename(test.name)}.comp.cmp.out"
        )
        with open(fpth, "w") as f:
            line = f"{'TOTIM':>15s}"
            line += f" {'CSUB':>15s}"
            line += f" {'MF':>15s}"
            line += f" {'DIFF':>15s}"
            f.write(line + "\n")
            for i in range(diff.shape[0]):
                line = f"{tc0['time'][i]:15g}"
                line += f" {tc['TCOMP3'][i]:15g}"
                line += f" {tc0['TCOMP3'][i]:15g}"
                line += f" {diff[i]:15g}"
                f.write(line + "\n")

        if diffmax > dtol:
            test.success = False
            msg += f"exceeds {dtol}"
            assert diffmax < dtol, msg
        else:
            test.success = True
            print("    " + msg)

    # compare budgets
    cbc_compare(test)


def cbc_compare(test):
    # open cbc file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.cbc")
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double")

    # build list of cbc data to retrieve
    avail = cobj.get_unique_record_names()
    cbc_bud = []
    bud_lst = []
    for t in avail:
        if isinstance(t, bytes):
            t = t.decode()
        t = t.strip()
        if paktest in t.lower():
            cbc_bud.append(t)
            bud_lst.append(f"{t}_IN")
            bud_lst.append(f"{t}_OUT")

    # get results from listing file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.lst")
    budl = flopy.utils.Mf6ListBudget(fpth)
    names = list(bud_lst)
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]
    d = np.recarray(nbud, dtype=dtype)
    for key in bud_lst:
        d[key] = 0.0

    # get data from cbc dile
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

    diff = np.zeros((nbud, len(bud_lst)), dtype=float)
    for i, key in enumerate(bud_lst):
        diff[:, i] = d0[key] - d[key]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-budget difference ({diffmax}) "

    # write summary
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.bud.cmp.out")
    with open(fpth, "w") as f:
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

    if diffmax > budtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)


@pytest.mark.slow
@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        htol=htol[idx],
        compare="mf6_regression",
        verbose=False,
    )
    test.run()

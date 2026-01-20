import os

import flopy
import numpy as np
import pytest
from conftest import try_get_target
from flopy.utils.compare import compare_heads
from framework import TestFramework

cases = ["csub_zdisp01"]
cmppth = "mfnwt"
htol = [None for _ in range(len(cases))]
dtol = 1e-3
budtol = 1e-2
bud_lst = [
    "STO-SS_IN",
    "STO-SS_OUT",
    "STO-SY_IN",
    "STO-SY_OUT",
    "CSUB-CGELASTIC_IN",
    "CSUB-CGELASTIC_OUT",
    "CSUB-ELASTIC_IN",
    "CSUB-ELASTIC_OUT",
    "CSUB-INELASTIC_IN",
    "CSUB-INELASTIC_OUT",
    "CSUB-WATERCOMP_IN",
    "CSUB-WATERCOMP_OUT",
]

# static model data
# temporal discretization
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
# tsmult = [1.0] + [1.0 for _ in range(nper - 1)]
steady = [True] + [False for _ in range(nper - 1)]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# spatial discretization data
nlay, nrow, ncol = 3, 20, 20
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1000.0, 2000.0
top = 0.0
botm = [-100, -150.0, -350.0]
zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
strt = 0.0
hnoflo = 1e30
hdry = -1e30

# create idomain and ibound
idomain = np.ones((nlay, nrow, ncol), dtype=np.int32)
idomain[0, 10:, :] = 0
idomain[1, 0:5, :] = 0
idomain[1, 15:, :] = 0
idomain[2, 0:10, :] = 0
iex = np.zeros((nlay, nrow, ncol), dtype=np.int32)
iex[idomain == 0] = 1

# calculate hk
hk1fact = 1.0 / 50.0
hk1 = 0.5 * hk1fact
# hk1[0, :] = 1000. * hk1fact
# hk1[-1, :] = 1000. * hk1fact
# hk1[:, 0] = 1000. * hk1fact
# hk1[:, -1] = 1000. * hk1fact
hk = [20.0, hk1, 5.0]

# calculate vka
vka = [1e6, 7.5e-5, 1e6]

# layer 1 is convertible
laytyp = [1, 0, 0]

# solver options
nouter, ninner = 500, 300
hclose, rclose, relax = 1e-9, 1e-6, 1.0
newtonoptions = "NEWTON"
imsla = "BICGSTAB"

# chd data
c = []
c6 = []
ccol = [j for j in range(ncol)]
for j in ccol:
    c.append([0, 0, j, strt, strt])
    c6.append([(0, 0, j), strt])
cd = {0: c}
cd6 = {0: c6}
maxchd = len(cd[0])

# drain data
dr = []
dr6 = []
drh = strt - 1.0
drc = 10.0
for j in ccol:
    dr.append([2, nrow - 1, j, drh, drc])
    dr6.append([(2, nrow - 1, j), drh, drc])
drd = {0: dr}
drd6 = {0: dr6}
maxdrd = len(drd[0])

# pumping well data
wrp = [12, 12, 13, 13]
wcp = [9, 10, 9, 10]
wq = [-14000.0, -8000.0, -5000.0, -3000.0]
d = []
d6 = []
for r, c, q in zip(wrp, wcp, wq):
    d.append([2, r, c, q])
    d6.append([(2, r, c), q])
wd = {1: d}
wd6 = {1: d6}
maxwel = len(wd[1])
maxwel = len(wd[1])

# storage and compaction data
# ske = [6e-4, 3e-4, 6e-4]
# ss = [3e-6, 0., 3e-6]
ss = [0.0, 0.0, 0.0]
void = 0.82
theta = void / (1.0 + void)

# static ibc and sub data
sgm = 0.0
sgs = 0.0
omega = 1.0

# no delay bed data
nndb = 3
lnd = [0, 1, 2]
hc = -7.0
thicknd0 = [15.0, 50.0, 30.0]
ccnd0 = [6e-4, 3e-4, 6e-4]
crnd0 = [6e-6, 3e-6, 6e-6]
sfv = []
sfe = []
for k in range(nlay):
    sfv.append(ccnd0[k] * thicknd0[k])
    sfe.append(crnd0[k] * thicknd0[k])

# ibc packagedata container counter
sub6 = []
ibcno = 0

# create no delay bed packagedata entries
if nndb > 0:
    cdelays = "nodelay"
    for kdx, k in enumerate(lnd):
        for i in range(nrow):
            for j in range(ncol):
                # skip constant head cells
                if idomain[k, i, j] == 0:
                    continue
                tag = f"{k + 1:02d}_{i + 1:02d}_{j + 1:02d}"
                # create nodelay entry
                # no delay beds
                b = thicknd0[kdx]
                d = [
                    ibcno,
                    (k, i, j),
                    cdelays,
                    hc,
                    b,
                    1.0,
                    ccnd0[kdx],
                    crnd0[kdx],
                    theta,
                    999.0,
                    -999.0,
                    tag,
                ]
                sub6.append(d)
                ibcno += 1

# create delay bed packagedata entries and coarse-grained materia storage
ske_scaled = []
# create S for aquifer and no-delay beds
for k in range(nlay):
    sst = (zthick[k] - thicknd0[k]) * ss[k] / zthick[k]
    ske_scaled.append(sst)

maxcsub = len(sub6)

# sub output data
ds15 = [0, 0, 0, 2052, 0, 0, 0, 2053, 0, 0, 0, 0]
ds16 = [0, nper - 1, 0, nstp[-1] - 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1]


# variant SUB package problem 3
def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
    elevs = [top] + botm

    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions=newtonoptions, save_flows=True
    )

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
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
        idomain=idomain,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf, save_flows=False, icelltype=laytyp, k=hk, k33=vka
    )
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0,
        sy=0,
        storagecoefficient=None,
        steady_state={0: True},
        transient={1: True},
    )

    # csub files
    opth = f"{name}.csub.obs"
    ibcsv = f"{name}.ib.strain.csv"
    skcsv = f"{name}.sk.strain.csv"
    copth = f"{name}.compaction.gridbin"
    zopth = f"{name}.zdisplacement.gridbin"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        boundnames=True,
        head_based=True,
        specified_initial_interbed_state=True,
        effective_stress_lag=True,
        save_flows=True,
        strainib_filerecord=ibcsv,
        straincg_filerecord=skcsv,
        compaction_filerecord=copth,
        zdisplacement_filerecord=zopth,
        ninterbeds=maxcsub,
        beta=0.0,
        cg_ske_cr=ss,
        packagedata=sub6,
    )
    orecarray = {}
    oloc = (2, wrp[0], wcp[0])
    ibloc = (449,)
    orecarray["csub_obs.csv"] = [
        ("tcomp3", "interbed-compaction", ibloc),
        ("sk-tcomp3", "coarse-compaction", oloc),
        ("ibi-tcomp3", "inelastic-compaction", ibloc),
        ("ibe-tcomp3", "elastic-compaction", ibloc),
    ]
    csub_obs_package = csub.obs.initialize(
        filename=opth, digits=10, print_input=True, continuous=orecarray
    )

    # drain
    drn = flopy.mf6.ModflowGwfdrn(gwf, maxbound=maxdrd, stress_period_data=drd6)

    # wel file
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=maxwel,
        stress_period_data=wd6,
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=maxchd, stress_period_data=cd6, save_flows=False
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

    # build MODFLOW-NWT files
    cpth = cmppth
    ws = os.path.join(test.workspace, cpth)
    mc = flopy.modflow.Modflow(
        name,
        model_ws=ws,
        version=cpth,
        exe_name=try_get_target(test.targets, "mfnwt"),
    )
    dis = flopy.modflow.ModflowDis(
        mc,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        nper=nper,
        perlen=perlen,
        nstp=nstp,
        tsmult=tsmult,
        steady=steady,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    bas = flopy.modflow.ModflowBas(
        mc, ibound=idomain, strt=strt, hnoflo=hnoflo, stoper=0.01
    )
    upw = flopy.modflow.ModflowUpw(
        mc,
        laytyp=laytyp,
        ipakcb=1001,
        hk=hk,
        vka=vka,
        ss=ske_scaled,
        sy=0.0,
        hdry=hdry,
    )
    sub = flopy.modflow.ModflowSub(
        mc,
        ndb=0,
        nndb=nndb,
        isuboc=1,
        ln=lnd,
        hc=hc,
        sfe=sfe,
        sfv=sfv,
        ids15=ds15,
        ids16=ds16,
    )
    chd = flopy.modflow.ModflowChd(mc, stress_period_data=cd)
    drn = flopy.modflow.ModflowDrn(mc, stress_period_data=drd)
    wel = flopy.modflow.ModflowWel(mc, stress_period_data=wd)
    oc = flopy.modflow.ModflowOc(
        mc,
        stress_period_data=None,
        save_every=1,
        save_types=["print head", "save head", "save budget"],
    )
    fluxtol = (float(nlay * nrow * ncol) - 4.0) * rclose
    nwt = flopy.modflow.ModflowNwt(
        mc,
        headtol=hclose,
        fluxtol=fluxtol,
        maxiterout=nouter,
        linmeth=2,
        maxitinner=ninner,
        unitnumber=132,
        options="SPECIFIED",
        backflag=0,
        idroptol=0,
    )

    return sim, mc


def check_output(idx, test):
    # MODFLOW 6 total compaction results
    fpth = os.path.join(test.workspace, "csub_obs.csv")
    try:
        tc = np.genfromtxt(fpth, names=True, delimiter=",")
    except:
        assert False, f'could not load data from "{fpth}"'

    # MODFLOW-2005 total compaction results
    fn = f"{os.path.basename(test.name)}.total_comp.hds"
    fpth = os.path.join(test.workspace, "mfnwt", fn)
    try:
        sobj = flopy.utils.HeadFile(fpth, text="LAYER COMPACTION", verbose=False)
        tc0 = sobj.get_ts((2, wrp[0], wcp[0]))
    except:
        assert False, f'could not load data from "{fpth}"'

    # calculate maximum absolute error
    diff = tc["TCOMP3"] - tc0[:, 1]
    diffmax = np.abs(diff).max()
    msg = f"maximum absolute total-compaction difference ({diffmax}) "

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
    else:
        test.success = True
        print("    " + msg)

    # get results from listing file
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.lst")
    budl = flopy.utils.Mf6ListBudget(fpth)
    names = list(bud_lst)
    d0 = budl.get_budget(names=names)[0]
    dtype = d0.dtype
    nbud = d0.shape[0]

    # get results from cbc file
    cbc_bud = [
        "STO-SS",
        "STO-SY",
        "CSUB-CGELASTIC",
        "CSUB-ELASTIC",
        "CSUB-INELASTIC",
        "CSUB-WATERCOMP",
    ]
    d = np.recarray(nbud, dtype=dtype)
    for key in bud_lst:
        d[key] = 0.0
    fpth = os.path.join(test.workspace, f"{os.path.basename(test.name)}.cbc")
    cobj = flopy.utils.CellBudgetFile(fpth, precision="double", verbose=False)
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

    # compare z-displacement data
    fpth1 = os.path.join(
        test.workspace,
        f"{os.path.basename(test.name)}.zdisplacement.gridbin",
    )
    fpth2 = os.path.join(test.workspace, cmppth, "csub_zdisp01.vert_disp.hds")
    text1 = "CSUB-ZDISPLACE"
    text2 = "Z DISPLACEMENT"
    fout = os.path.join(
        test.workspace,
        f"{os.path.basename(test.name)}.z-displacement.bin.out",
    )
    success_tst = compare_heads(
        None,
        None,
        text=text1,
        text2=text2,
        outfile=fout,
        files1=fpth1,
        files2=fpth2,
        difftol=True,
        verbose=True,
        exarr=iex,
    )
    msg = f"z-displacement comparison success = {success_tst}"
    if success_tst:
        test.success = True
        print(msg)
    else:
        test.success = False
        assert success_tst, msg


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
    )
    test.run()

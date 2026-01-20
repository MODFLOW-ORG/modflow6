import flopy
import pytest
from framework import TestFramework

# fmt: off
obs_names = [
    "delay-flowtop", "delay-flowtop",
    "delay-flowbot", "delay-flowbot",
    "delay-head", "delay-head",
    "delay-gstress", "delay-gstress", 
    "delay-estress", "delay-estress", 
    "delay-preconstress", "delay-preconstress", 
    "delay-compaction", "delay-compaction",
    "delay-thickness", "delay-thickness",
    "delay-theta", "delay-theta",
    "csub", "csub",
    "inelastic-csub", "inelastic-csub", 
    "elastic-csub", "elastic-csub",
    "interbed-compaction", "interbed-compaction",
    "inelastic-compaction", "inelastic-compaction",
    "elastic-compaction", "elastic-compaction",
    ]
# fmt: off
boundname = [
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    ]
# fmt: off
test_fail = [
    False, False,
    False, False,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, True,
    False, False,
    False, False,
    False, False,
    False, True,
    False, True,
    False, True,
    ]
cases = [f"csub_obs{idx + 1:02d}" for idx, _ in enumerate(obs_names)]

paktest = "csub"
budtol = 1e-2
ndcell = [19] * len(cases)

# static model data
# spatial discretization
nlay, nrow, ncol = 1, 1, 3
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1.0, 1.0
top = 0.0
botm = [-100.0]

# temporal discretization
nper = 1
perlen = [1000.0 for _ in range(nper)]
nstp = [100 for _ in range(nper)]
tsmult = [1.05 for _ in range(nper)]
steady = [False for _ in range(nper)]

strt = 0.0
strt6 = 1.0
hnoflo = 1e30
hdry = -1e30
hk = 1e6
laytyp = [0]
S = 1e-4
sy = 0.0

nouter, ninner = 1000, 300
hclose, rclose, relax = 1e-6, 1e-6, 0.97

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

ib = 1

c = []
c6 = []
for j in range(0, ncol, 2):
    c.append([0, 0, j, strt, strt])
    c6.append([(0, 0, j), strt])
cd = {0: c}
cd6 = {0: c6}

# sub data
ndb = 1
nndb = 0
cc = 100.0
cr = 1.0
void = 0.82
theta = void / (1.0 + void)
kv = 0.025
sgm = 0.0
sgs = 0.0
ini_stress = 1.0
thick = [1.0]
sfe = cr * thick[0]
sfv = cc * thick[0]
lnd = [0]
ldnd = [0]
dp = [[kv, cr, cc]]
ss = S / (100.0 - thick[0])

ds15 = [0, 0, 0, 2052, 0, 0, 0, 0, 0, 0, 0, 0]
ds16 = [0, 0, 0, 100, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]


def get_model(idx, ws):
    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution
    ims = flopy.mf6.ModflowIms(
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
    )

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name)

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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=hk)
    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        storagecoefficient=True,
        transient={0: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf, maxbound=len(c6), stress_period_data=cd6, save_flows=False
    )

    # csub files
    sub6 = [
        [
            0,
            (0, 0, 1),
            "delay",
            ini_stress,
            thick[0],
            1.0,
            230.258658761733000,
            2.302586587617330,
            theta,
            kv,
            ini_stress,
        ]
    ]
    bname = "interbed"
    if boundname[idx]:
        sub6[0].append(bname)
        obs_idx = ("obs_value", obs_names[idx], bname)
    else:
        obs_idx = ("obs_value", obs_names[idx], (0,), (0,))

    opth = f"{name}.csub.obs"
    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        boundnames=boundname[idx],
        head_based=False,
        print_input=True,
        save_flows=True,
        ndelaycells=ndcell[idx],
        ninterbeds=1,
        beta=0.0,
        cg_ske_cr=ss,
        packagedata=sub6,
    )
    orecarray = {}
    orecarray["csub_obs.csv"] = [obs_idx]
    csub_obs_package = csub.obs.initialize(
        filename=opth, digits=10, print_input=True, continuous=orecarray
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        printrecord=[("BUDGET", "ALL")],
    )

    return sim


def build_models(idx, test):
    sim = get_model(idx, test.workspace)

    return sim, None


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
        xfail=test_fail[idx],
    )
    test.run()

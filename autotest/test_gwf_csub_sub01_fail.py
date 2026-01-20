#
# Modified version of sub example 1 where the interbed thickness exceeds the
# layer thickness and causes the model to terminate with an error
#


import flopy
import pytest
from framework import TestFramework

cases = ["csub_sub01"]
paktest = "csub"
ndcell = [19] * len(cases)

# static model data
# spatial discretization
nlay, nrow, ncol = 1, 1, 3
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1.0, 1.0
top = 0.0
botm = [-10.0]

# temporal discretization
nper = 1
perlen = [1000.0] * nper
nstp = [100] * nper
tsmult = [1.05] * nper
steady = [False] * nper

strt = 0.0
strt6 = 1.0
hk = 1e6
laytyp = [0]
S = 1e-4
sy = 0.0

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
ninterbeds = 3
cc = 100.0
cr = 1.0
void = 0.82
theta = void / (1.0 + void)
kv = 0.025
sgm = 0.0
sgs = 0.0
ini_stress = 1.0
thick = [1.1, 2.2, 5.5]
rnb = [3.3, 2.1, 1.0]
cdelay = ["delay", "delay", "nodelay"]
boundnames = [f"interbed_{i + 1}" for i in range(ninterbeds)]
sfe = [cr * b for b in thick]
sfv = [cc * b for b in thick]


def get_model(idx, ws):
    name = cases[idx]

    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        version="mf6",
        exe_name="mf6",
        sim_ws=ws,
        print_input=True,
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create iterative model solution
    ims = flopy.mf6.ModflowIms(sim, complexity="simple")

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
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=laytyp, k=hk, k33=hk)

    # storage
    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=laytyp,
        ss=0.0,
        sy=sy,
        storagecoefficient=True,
        transient={0: True},
    )

    # chd files
    chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(
        gwf,
        maxbound=len(c6),
        stress_period_data=cd6,
    )

    # csub files
    sub6 = []
    for n in range(ninterbeds):
        sub6.append(
            (
                n,
                (0, 0, 1),
                cdelay[n],
                ini_stress,
                thick[n],
                rnb[n],
                sfv[n],
                sfe[n],
                theta,
                kv,
                ini_stress,
                boundnames[n],
            )
        )

    csub = flopy.mf6.ModflowGwfcsub(
        gwf,
        print_input=True,
        head_based=True,
        boundnames=True,
        save_flows=True,
        effective_stress_lag=True,
        ndelaycells=ndcell[idx],
        ninterbeds=ninterbeds,
        beta=0.0,
        cg_ske_cr={"iprn": 1, "data": 1e-5},
        packagedata=sub6,
    )

    return sim


def build_model(idx, test):
    # build MODFLOW 6 files
    sim = get_model(idx, test.workspace)

    return sim, None


def check_output(idx, test):
    buff = test.buffs[0]
    assert any("Coarse grained material thickness is less than zero" in l for l in buff)


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_model(idx, t),
        targets=targets,
        check=lambda t: check_output(idx, t),
        xfail=True,
    )
    test.run()

import pathlib as pl

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["csub_zdisp02"]
htol = [None for _ in range(len(cases))]
dtol = 1e-3
budtol = 1e-2

# static model data
# temporal discretization
nper = 31
perlen = [1.0] + [365.2500000 for _ in range(nper - 1)]
nstp = [1] + [6 for _ in range(nper - 1)]
tsmult = [1.0] + [1.3 for _ in range(nper - 1)]
tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))

# spatial discretization data
nlay, nrow, ncol = 3, 6, 6
shape3d = (nlay, nrow, ncol)
size3d = nlay * nrow * ncol
delr, delc = 1000.0, 1000.0
top = 0.0
botm = [-100, -150.0, -350.0]
zthick = [top - botm[0], botm[0] - botm[1], botm[1] - botm[2]]
strt = 0.0

# create idomain and ibound
idomain = np.ones(shape3d, dtype=np.int32)
idomain[0, :, :] = np.array(
    [
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1, 1, 1],
        [0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 0, 0],
        [1, 1, 1, 0, 0, 0],
    ]
)
idomain[1, :, :] = np.flip(idomain[0], axis=1)
idomain[2] = idomain[0] * idomain[1]

# calculate hk
hk1fact = 1.0 / 50.0
hk1 = 0.5 * hk1fact
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


# pumping well data
spd_wel = [(2, 3, 2, -1400.0)]
maxwel = len(spd_wel)

# storage and compaction data
ss = [0.0, 0.0, 0.0]
void = 0.82
theta = void / (1.0 + void)

# static ibc and sub data
sgm = 0.0
sgs = 0.0
omega = 1.0

# no delay bed data
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

# no delay bed packagedata entries
sub6 = []
cdelays = "nodelay"
ibcno = 0
for kdx, k in enumerate(lnd):
    for i in range(nrow):
        for j in range(ncol):
            # skip inactive cells
            if idomain[k, i, j] == 0:
                continue
            tag = f"{k + 1:02d}_{i + 1:02d}_{j + 1:02d}"
            # create nodelay entry
            b = thicknd0[kdx]
            d = (
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
            )
            sub6.append(d)
            ibcno += 1

# create coarse-grained material storage
ske_scaled = []
for k in range(nlay):
    sst = (zthick[k] - thicknd0[k]) * ss[k] / zthick[k]
    ske_scaled.append(sst)


def build_models(idx, test):
    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )

    flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    flopy.mf6.ModflowIms(
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

    # create gwf model
    gwf = flopy.mf6.ModflowGwf(
        sim, modelname=name, newtonoptions=newtonoptions, save_flows=True
    )

    flopy.mf6.ModflowGwfdis(
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

    flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{name}.ic")

    flopy.mf6.ModflowGwfnpf(gwf, save_flows=False, icelltype=laytyp, k=hk, k33=vka)

    flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=False,
        iconvert=laytyp,
        ss=0,
        sy=0,
        storagecoefficient=None,
        steady_state={0: True},
        transient={1: True},
    )

    copth = f"{name}.compaction.gridbin"
    zopth = f"{name}.zdisplacement.gridbin"
    cvgpth = f"{name}.csub.convergence.csv"
    flopy.mf6.ModflowGwfcsub(
        gwf,
        boundnames=True,
        head_based=True,
        update_material_properties=True,
        specified_initial_interbed_state=True,
        effective_stress_lag=False,
        save_flows=True,
        compaction_filerecord=copth,
        zdisplacement_filerecord=zopth,
        package_convergence_filerecord=cvgpth,
        ninterbeds=len(sub6),
        beta=0.0,
        cg_ske_cr=ss,
        packagedata=sub6,
    )

    flopy.mf6.ModflowGwfwel(
        gwf,
        print_input=True,
        print_flows=True,
        maxbound=maxwel,
        stress_period_data={1: spd_wel},
    )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "LAST"), ("BUDGET", "ALL")],
    )

    return sim, None


def check_output(idx, test):
    with open(pl.Path(test.workspace) / "mfsim.lst", "r") as f:
        lines = f.readlines()
    tag = "1. Package convergence data is requested but delay interbeds are not"
    found = False
    for line in lines[::-1]:
        if tag in line:
            found = True
            break
    assert found, f"WARNING: '{tag}' not found in mfsim.lst"

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwf = sim.get_model()
    idomain = gwf.dis.idomain.array

    times = gwf.csub.output.compaction().get_times()
    totim = times[-1]

    compaction = gwf.csub.output.compaction().get_data(totim=totim)
    compaction[compaction == 1e30] = 0.0
    zdis_calc = np.zeros(compaction.shape, dtype=float)
    zdis_calc[2] = compaction[2]
    for k in (1, 0):
        zdis_calc[k] = zdis_calc[k + 1] + compaction[k]
    zdis_calc[idomain < 1] = 0.0

    zpth = pl.Path(test.workspace) / "csub_zdisp02.zdisplacement.gridbin"
    zobj = flopy.utils.HeadFile(zpth, text="CSUB-ZDISPLACE")
    zdis = zobj.get_data(totim=totim)
    zdis[zdis == 1e30] = 0.0

    assert np.allclose(zdis, zdis_calc), (
        "Calculated z-displacement is not equal to simulated z-displacement"
    )


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

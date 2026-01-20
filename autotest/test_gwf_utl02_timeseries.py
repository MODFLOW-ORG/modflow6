import flopy
import pytest
from framework import TestFramework

cases = ["ts01"]


def build_models(idx, test):
    nlay, nrow, ncol = 1, 3, 3
    nper = 2
    perlen = [1.0, 14966]
    nstp = 1
    tsmult = 1.0
    steady = [True, False]
    lenx = 300.0
    delr = delc = lenx / float(nrow)
    botm = -1.0
    hk = 1.0

    nouter, ninner = 100, 300
    hclose, rclose, relax = 1e-6, 1e-3, 1.0

    tdis_rc = []
    for i in range(nper):
        tdis_rc.append((perlen[i], nstp, tsmult))

    name = cases[idx]

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)

    # create gwf model
    gwf = flopy.mf6.MFModel(
        sim,
        model_type="gwf6",
        modelname=name,
        model_nam_file=f"{name}.nam",
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
        linear_acceleration="CG",
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
        top=0.0,
        botm=botm,
        idomain=1,
        filename=f"{name}.dis",
    )

    # initial conditions
    ic = flopy.mf6.ModflowGwfic(gwf, strt=0.0, filename=f"{name}.ic")

    # node property flow
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        icelltype=0,
        k=hk,
        k33=hk,
        filename=f"{name}.npf",
    )

    # chd files
    chdlist0 = []
    chdlist0.append([(0, 0, 0), 1.0])
    chdlist0.append([(nlay - 1, nrow - 1, ncol - 1), 0.0])

    chdspdict = {0: chdlist0}
    chd = flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chdspdict,
        save_flows=False,
        filename=f"{name}.chd",
    )

    # wel files
    wellist1 = []
    wellist1.append([(0, 2, 2), "ts01"])
    wellist1.append([(0, 2, 2), "ts02"])
    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        pname="wel",
        print_input=True,
        print_flows=True,
        stress_period_data={1: wellist1},
    )
    # ts_filerecord='well-rates.ts')

    # well ts package
    ts_recarray = [
        (0.0, 0.0, 0.0),
        (700.0, 3.0e30, -0.5),
        (31775.0, -0.5, 3.0e30),
        (41272.0, 3.0e30, 3.0e30),
        (15000.0, 0.0, 0.0),
    ]

    filename = name + ".wel.ts"
    time_series_namerecord = [("ts01", "ts02")]
    interpolation_methodrecord = [("linear", "linear")]
    wel.ts.initialize(
        filename=filename,
        timeseries=ts_recarray,
        time_series_namerecord=time_series_namerecord,
        interpolation_methodrecord=interpolation_methodrecord,
    )

    # output control
    oc = flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        filename=f"{name}.oc",
    )

    return sim, None


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        targets=targets,
    )
    test.run()

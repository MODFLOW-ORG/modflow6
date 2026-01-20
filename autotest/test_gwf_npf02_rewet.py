import os

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["npf02_hreweta", "npf02_hrewetb", "npf02_hrewetc", "npf02_hrewetd"]
ncols = [[15], [10, 5], [15], [10, 5]]
nlays = [1, 1, 3, 3]

# static model data
nrow = 10
nper = 2
perlen = [1.0, 1.0]
nstp = [1, 1]
tsmult = [1.0, 1.0]

lenx = 15.0 * 500.0
leny = 10.0 * 500.0
delr = lenx / float(ncols[0][0])
delc = leny / float(nrow)
top = 150.0
strt = -40.0
hk = 10.0

# rewetting
rewet_record = [("WETFCT", 1.0, "IWETIT", 1, "IHDWET", 1)]
wetdry = -0.001  # [0.001, 0.001, 0.001]

# left chd boundary for each stress period
hbndl = [100.0, 25]

nouter, ninner = 1000, 100
hclose, rclose, relax = 1e-1, 0.01, 1.0

tdis_rc = []
for i in range(nper):
    tdis_rc.append((perlen[i], nstp[i], tsmult[i]))


def get_local_data(idx):
    ncolst = ncols[idx]
    nmodels = len(ncolst)
    mnames = []
    for jdx in range(nmodels):
        mname = f"gwf{jdx}"
        mnames.append(mname)
    return ncolst, nmodels, mnames


def build_models(idx, test):
    name = cases[idx]
    nlay = nlays[idx]

    if nlay == 1:
        botm = [-50.0]
    elif nlay == 3:
        botm = [50.0, 0.0, -50.0]

    c6left = []
    c6right = []
    vl = hbndl[0]
    vr = strt
    for k in range(nlay):
        for i in range(nrow):
            if botm[k] < vl:
                c6left.append([(k, i, 0), vl])
            if botm[k] < vr:
                c6right.append([(k, i, ncols[idx][-1] - 1), vr])
    cd6left = {0: c6left}
    cd6right = {0: c6right}
    c6left = []
    vl = hbndl[1]
    for k in range(nlay):
        for i in range(nrow):
            if botm[k] < vl:
                c6left.append([(k, i, 0), vl])
    cd6left[1] = c6left

    # build MODFLOW 6 files
    ws = test.workspace
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    # create tdis package
    tdis = flopy.mf6.ModflowTdis(sim, time_units="DAYS", nper=nper, perioddata=tdis_rc)
    # set ims csv files
    csv0 = f"{name}.outer.ims.csv"
    csv1 = f"{name}.inner.ims.csv"

    # create iterative model solution and register the gwf model with it
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        csv_outer_output_filerecord=csv0,
        csv_inner_output_filerecord=csv1,
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

    # set local data for this model
    ncolst, nmodels, mnames = get_local_data(idx)

    sim.register_ims_package(ims, mnames)

    # create exchange file, if needed
    if nmodels > 1:
        exchd = []
        for k in range(nlay):
            for i in range(nrow):
                t = []
                t.append((k, i, ncolst[0] - 1))
                t.append((k, i, 0))
                t.append(1)
                t.append(delr / 2.0)
                t.append(delr / 2.0)
                t.append(delc)
                exchd.append(t)
        excf = flopy.mf6.ModflowGwfgwf(
            sim,
            exgtype="GWF6-GWF6",
            nexg=len(exchd),
            exgmnamea=mnames[0],
            exgmnameb=mnames[1],
            exchangedata=exchd,
        )

    # create gwf model
    for jdx in range(nmodels):
        mname = mnames[jdx]

        gwf = flopy.mf6.ModflowGwf(sim, modelname=mname, model_nam_file=f"{mname}.nam")

        dis = flopy.mf6.ModflowGwfdis(
            gwf,
            nlay=nlay,
            nrow=nrow,
            ncol=ncolst[jdx],
            delr=delr,
            delc=delc,
            top=top,
            botm=botm,
            filename=f"{mname}.dis",
        )

        # initial conditions
        ic = flopy.mf6.ModflowGwfic(gwf, strt=strt, filename=f"{mname}.ic")

        # node property flow
        npf = flopy.mf6.ModflowGwfnpf(
            gwf,
            save_flows=False,
            rewet_record=rewet_record,
            icelltype=1,
            k=hk,
            wetdry=wetdry,
        )

        # chd files
        if jdx == 0:
            fn = f"{mname}.chd1.chd"
            chd1 = flopy.mf6.modflow.ModflowGwfchd(
                gwf,
                stress_period_data=cd6left,
                save_flows=False,
                filename=fn,
                pname="chd1",
                print_input=True,
            )
        if jdx == nmodels - 1:
            fn = f"{mname}.chd2.chd"
            chd2 = flopy.mf6.modflow.ModflowGwfchd(
                gwf,
                stress_period_data=cd6right,
                save_flows=False,
                filename=fn,
                pname="chd2",
                print_input=True,
            )

        # output control
        oc = flopy.mf6.ModflowGwfoc(
            gwf,
            budget_filerecord=f"{mname}.cbc",
            head_filerecord=f"{mname}.hds",
            headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
            saverecord=[("HEAD", "LAST")],
            printrecord=[("HEAD", "LAST"), ("BUDGET", "LAST")],
        )

    return sim, None


def check_output(idx, test):
    hdata01lay = [
        [
            1.000000000000000000e02,
            9.491194679807708212e01,
            8.963852725425633139e01,
            8.415783778939784554e01,
            7.844327180838388358e01,
            7.246196989134719502e01,
            6.617253575516674857e01,
            5.952154961171697778e01,
            5.243800230005604845e01,
            4.482387907284233819e01,
            3.653700567841701030e01,
            2.735652405614735727e01,
            1.690257501779389671e01,
            4.399393234039389533e00,
            -4.000000000000000000e01,
        ],
        [
            2.500000000000000000e01,
            2.236586340188107513e01,
            1.963190756499453116e01,
            1.678583852445553148e01,
            1.381261073936328998e01,
            1.069347358389677538e01,
            7.404548318470112633e00,
            3.914611143252753056e00,
            1.814550305069945468e-01,
            -3.854475465269350920e00,
            -8.282216319385574010e00,
            -1.324637250326845006e01,
            -1.901458294740395516e01,
            -2.622133464488191024e01,
            -4.000000000000000000e01,
        ],
    ]
    hdata03lay = [
        [
            1.000000000000000000e02,
            9.496635368428880497e01,
            8.974621521357816789e01,
            8.432145824563308167e01,
            7.866533338822696919e01,
            7.274409668492738490e01,
            6.651345136714729733e01,
            5.991038738691469945e01,
            5.278755954306133447e01,
            4.516004619735397796e01,
            3.696996838160407606e01,
            2.789393025475152399e01,
            1.752424435308674333e01,
            4.558574627023233461e00,
            -4.000000000000000000e01,
        ],
        [
            2.500000000000000000e01,
            2.237276762221516435e01,
            1.964557102731594540e01,
            1.680659328791805862e01,
            1.384042482892330028e01,
            1.072732402849822897e01,
            7.440599125420748194e00,
            3.939435701073743079e00,
            9.093214618820914807e-02,
            -3.940381765174159057e00,
            -8.354930726606033531e00,
            -1.330380168293016219e01,
            -1.905635937367108212e01,
            -2.625249959368398223e01,
            -4.000000000000000000e01,
        ],
    ]

    ncolst, nmodels, mnames = get_local_data(idx)
    nlay = nlays[idx]

    # make single head array
    ncolt = 0
    for ncol in ncolst:
        ncolt += ncol
    hval = np.zeros((nper, ncolt), dtype=float)
    imid = int(nrow / 2)

    for j in range(nmodels):
        fn = os.path.join(test.workspace, f"{mnames[j]}.hds")
        hobj = flopy.utils.HeadFile(fn)
        times = hobj.get_times()
        ioff = 0
        if j > 0:
            ioff += ncolst[j - 1]
        for n, t in enumerate(times):
            h = hobj.get_data(totim=t)[:, imid, :]
            h = h.reshape((nlay, ncolst[j]))
            ht = np.ones(h.shape[1], dtype=float) * 1e30
            for k in range(nlay):
                v = h[k, :].copy()
                kdx = (ht == 1e30) & (v != -1e30)
                ht[kdx] = v[kdx]
            i1 = ioff + ht.shape[0]
            hval[n, ioff:i1] = ht.copy()

    # # save results if the know results change slightly
    # fpth = os.path.join(sim.workspace, "results.dat")
    # np.savetxt(fpth, hval, delimiter=",")

    # known results
    if idx < 2:
        h0 = np.array(hdata01lay)
    else:
        h0 = np.array(hdata03lay)

    # calculate maximum absolute error
    diff = hval - h0
    diffmax = np.abs(diff).max()
    dtol = 1e-9
    msg = f"maximum absolute maw head difference ({diffmax}) "

    if diffmax > dtol:
        test.success = False
        msg += f"exceeds {dtol}"
        assert diffmax < dtol, msg
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
    )
    test.run()

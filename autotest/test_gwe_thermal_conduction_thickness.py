"""GWE LKE/SFE/MWE PACKAGEDATA thermal-conduction thickness (RBTHCND/FTHK)
must be > 0. Thickness is used as a bare divisor (ctherm = ktf * wa / s) in
the matrix formulation and budget-term routines, so a thickness <= 0 must
be rejected at PACKAGEDATA-read time to prevent a divide-by-zero.
"""

import flopy
import pytest


def run_mf6(argv, ws):
    import subprocess

    buff = []
    proc = subprocess.Popen(
        argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ws
    )
    result, _ = proc.communicate()
    if result is not None:
        c = result.decode("utf-8").rstrip("\r\n")
        print(f"{c}")
        buff.append(c)
    return proc.returncode, buff


def run_mf6_error(ws, exe, err_str_list):
    returncode, buff = run_mf6([exe], ws)
    msg = "mf terminated with error"
    if returncode != 0:
        if not isinstance(err_str_list, list):
            err_str_list = [err_str_list]
        for err_str in err_str_list:
            if any(err_str in s for s in buff):
                raise RuntimeError(msg)
            else:
                msg += " but did not print correct error message."
                msg += f'  Correct message should have been "{err_str}"'
                raise ValueError(msg)
    else:
        raise ValueError("mf6 terminated successfully but was expected to fail")


EXPECTED = "Specified thickness used for thermal conduction MUST BE >"


def test_lke_zero_thickness(function_tmpdir, targets):
    mf6 = targets["mf6"]
    name = "lke_zthk"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam")
    flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", complexity="SIMPLE", filename=f"{gwfname}.ims"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=5, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGwfic(gwf, strt=0.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=20.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), -0.5, 0.0], [(0, 0, 4), -0.5, 0.0]],
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )
    connlen = connwidth = 0.5
    flopy.mf6.ModflowGwflak(
        gwf,
        nlakes=1,
        noutlets=0,
        ntables=0,
        packagedata=[(0, -0.4, 3, 0.0)],
        connectiondata=[
            (0, 0, (0, 0, 1), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 1, (0, 0, 3), "HORIZONTAL", 0.0, 10, 10, connlen, connwidth),
            (0, 2, (0, 0, 2), "VERTICAL", 0.0, 10, 10, connlen, connwidth),
        ],
        perioddata=[(0, "STATUS", "CONSTANT"), (0, "STAGE", -0.4)],
        pname="LAK-1",
        auxiliary=["CONCENTRATION"],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{gwfname}.cbc")

    gwename = "gwe_" + name
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwename}.ims"), [gwe.name])
    flopy.mf6.ModflowGwedis(
        gwe, nlay=1, nrow=1, ncol=5, delr=1.0, delc=1.0, top=0.0, botm=-1.0
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=4183.0,
        density_water=999.728,
        latent_heat_vaporization=2500.0,
        heat_capacity_solid=800.0,
        density_solid=2650.0,
    )
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.2700)
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGweoc(gwe, budget_filerecord=f"{gwename}.cbc")
    # ktf, rbthcnd: rbthcnd=0.0 must be rejected
    flopy.mf6.ModflowGwelke(
        gwe,
        packagedata=[(0, 35.0, 0.5, 0.0, "mylake")],
        boundnames=True,
        lakeperioddata=[(0, "STATUS", "CONSTANT"), (0, "TEMPERATURE", 100.0)],
        flow_package_name="LAK-1",
        pname="LKE-1",
    )
    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(str(function_tmpdir), mf6, EXPECTED)


def test_sfe_zero_thickness(function_tmpdir, targets):
    mf6 = targets["mf6"]
    name = "sfe_zthk"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam")
    flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", complexity="SIMPLE", filename=f"{gwfname}.ims"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 5.0, 0.0], [(0, 0, 2), 4.0, 0.0]],
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )
    flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=1,
        packagedata=[
            [0, (0, 0, 1), 100.0, 5.0, 1e-3, 4.0, 1.0, 1e-5, 0.04, 0, 1.0, 0, 0.0]
        ],
        connectiondata=[[0]],
        perioddata={0: [[0, "status", "active"], [0, "inflow", 1.0]]},
        pname="SFR-1",
        auxiliary=["CONCENTRATION"],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{gwfname}.cbc")

    gwename = "gwe_" + name
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwename}.ims"), [gwe.name])
    flopy.mf6.ModflowGwedis(
        gwe, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[0.0]
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=4183.0,
        density_water=999.728,
        latent_heat_vaporization=2500.0,
        heat_capacity_solid=800.0,
        density_solid=2650.0,
    )
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.2700)
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGweoc(gwe, budget_filerecord=f"{gwename}.cbc")
    # ktf, rbthcnd: rbthcnd=0.0 must be rejected
    flopy.mf6.ModflowGwesfe(
        gwe,
        packagedata=[(0, 0.0, 0.5, 0.0)],
        reachperioddata={0: [(0, "STATUS", "ACTIVE")]},
        flow_package_name="SFR-1",
        pname="SFE-1",
    )
    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(str(function_tmpdir), mf6, EXPECTED)


def test_mwe_zero_thickness(function_tmpdir, targets):
    mf6 = targets["mf6"]
    name = "mwe_zthk"
    sim = flopy.mf6.MFSimulation(
        sim_name=name, exe_name=mf6, sim_ws=str(function_tmpdir)
    )
    flopy.mf6.ModflowTdis(sim, nper=1, perioddata=[(1.0, 1, 1.0)])
    gwfname = "gwf_" + name
    gwf = flopy.mf6.ModflowGwf(sim, modelname=gwfname, model_nam_file=f"{gwfname}.nam")
    flopy.mf6.ModflowIms(
        sim, print_option="SUMMARY", complexity="SIMPLE", filename=f"{gwfname}.ims"
    )
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[-100.0]
    )
    flopy.mf6.ModflowGwfnpf(gwf, k=1.0)
    flopy.mf6.ModflowGwfic(gwf, strt=5.0)
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=[[(0, 0, 0), 5.0, 0.0], [(0, 0, 2), 4.0, 0.0]],
        pname="CHD-1",
        auxiliary="CONCENTRATION",
    )
    flopy.mf6.ModflowGwfmaw(
        gwf,
        nmawwells=1,
        packagedata=[(0, 0.15, -100.0, 5.0, "thiem", 1, 0.0)],
        connectiondata=[(0, 0, (0, 0, 1), 10.0, -100.0, 1.0, 0.25)],
        perioddata={0: [(0, "status", "active"), (0, "rate", -1.0)]},
        pname="MAW-1",
        auxiliary=["CONCENTRATION"],
    )
    flopy.mf6.ModflowGwfoc(gwf, budget_filerecord=f"{gwfname}.cbc")

    gwename = "gwe_" + name
    gwe = flopy.mf6.ModflowGwe(sim, modelname=gwename, model_nam_file=f"{gwename}.nam")
    flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="SIMPLE",
        linear_acceleration="BICGSTAB",
        filename=f"{gwename}.ims",
    )
    sim.register_ims_package(sim.get_package(f"{gwename}.ims"), [gwe.name])
    flopy.mf6.ModflowGwedis(
        gwe, nlay=1, nrow=1, ncol=3, delr=100.0, delc=100.0, top=10.0, botm=[-100.0]
    )
    flopy.mf6.ModflowGweic(gwe, strt=0.0)
    flopy.mf6.ModflowGweadv(gwe, scheme="UPSTREAM")
    flopy.mf6.ModflowGweest(
        gwe,
        porosity=0.30,
        heat_capacity_water=4183.0,
        density_water=999.728,
        latent_heat_vaporization=2500.0,
        heat_capacity_solid=800.0,
        density_solid=2650.0,
    )
    flopy.mf6.ModflowGwecnd(gwe, xt3d_off=True, ktw=0.5918, kts=0.2700)
    flopy.mf6.ModflowGwessm(gwe, sources=[("CHD-1", "AUX", "CONCENTRATION")])
    flopy.mf6.ModflowGweoc(gwe, budget_filerecord=f"{gwename}.cbc")
    # ktf, fthk: fthk=0.0 must be rejected
    flopy.mf6.ModflowGwemwe(
        gwe,
        packagedata=[(0, 0.0, 0.5, 0.0, "well1")],
        boundnames=True,
        mweperioddata=[(0, "STATUS", "ACTIVE")],
        flow_package_name="MAW-1",
        pname="MWE-1",
    )
    flopy.mf6.ModflowGwfgwe(
        sim, exgtype="GWF6-GWE6", exgmnamea=gwfname, exgmnameb=gwename
    )
    sim.write_simulation()

    with pytest.raises(RuntimeError):
        run_mf6_error(str(function_tmpdir), mf6, EXPECTED)

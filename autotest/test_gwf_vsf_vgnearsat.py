"""
Test the near-saturation options of the Van Genuchten soil model in VSF:
KR_SMOOTHING (smooth relative permeability below saturation) and HAEVGN
(air-entry pressure head, Vogel et al., 2001).

A steady downward flux q passes through a tall silt loam column with the
water table at the bottom. Far above the water table the gradient is one,
so kr(psi) = q / Ks there. The simulated pressure head is compared with
the head that solves this equation for each variant of the kr curve, and,
for the air-entry variant, the saturation with the rescaled retention
curve.
"""

import flopy
import numpy as np
import pytest
from framework import TestFramework

cases = ["vgns-plain", "vgns-smooth", "vgns-airentry"]
kr_smoothing = [None, 2.0, None]
haevgn = [None, None, -2.0]

# silt loam (Carsel and Parrish, 1988), cm and s
ks, alpha, nvg = 1.25e-4, 0.020, 1.41
porosity, satres = 0.45, 0.067 / 0.45
kr_target = 0.7
q = kr_target * ks

nlay = 100
delz = 2.0
top = nlay * delz


def se_vg(psi):
    m = 1.0 - 1.0 / nvg
    return np.where(psi >= 0.0, 1.0, (1.0 + (alpha * np.abs(psi)) ** nvg) ** -m)


def kr_vg(psi):
    m = 1.0 - 1.0 / nvg
    se = se_vg(psi)
    return np.sqrt(se) * (1.0 - (1.0 - se ** (1.0 / m)) ** m) ** 2


def kr_curve(psi, idx):
    if haevgn[idx] is not None:
        hae = haevgn[idx]
        return 1.0 if psi >= hae else kr_vg(psi) / kr_vg(hae)
    if kr_smoothing[idx] is not None and psi > -kr_smoothing[idx]:
        d = kr_smoothing[idx]
        if psi >= 0.0:
            return 1.0
        eps = 1.0e-7
        k0 = kr_vg(-d)
        m0 = (kr_vg(-d + eps) - kr_vg(-d - eps)) / (2.0 * eps) * d
        t = (psi + d) / d
        return (
            (2 * t**3 - 3 * t**2 + 1) * k0
            + (t**3 - 2 * t**2 + t) * m0
            + (-2 * t**3 + 3 * t**2)
        )
    return kr_vg(psi)


def psi_unit_gradient(idx):
    """Pressure head where kr = q / Ks, by bisection (kr increases with psi)"""
    lo, hi = -1000.0, 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if kr_curve(mid, idx) < kr_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build_models(idx, test):
    name = cases[idx]
    botm = [top - (k + 1) * delz for k in range(nlay)]
    zc = np.array(botm) + 0.5 * delz

    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=test.workspace, exe_name="mf6")
    flopy.mf6.ModflowTdis(sim, time_units="SECONDS", perioddata=[(1.0, 1, 1.0)])
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1.0e-9,
        outer_maximum=500,
        inner_maximum=300,
        inner_dvclose=1.0e-10,
        rcloserecord=1.0e-12,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=1, ncol=1, delr=10.0, delc=10.0, top=top, botm=botm
    )
    # start with a pressure head of -10 cm everywhere
    flopy.mf6.ModflowGwfic(gwf, strt=zc - 10.0)
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=ks, k33=ks)
    flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    opts = {}
    if kr_smoothing[idx] is not None:
        opts["kr_smoothing"] = kr_smoothing[idx]
    if haevgn[idx] is not None:
        opts["haevgn"] = haevgn[idx]
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        soil_model="vangenuchten",
        porosity=porosity,
        satres=satres,
        alphavgn=alpha,
        nvgn=nvg,
        sat_filerecord=f"{name}.sat",
        phead_filerecord=f"{name}.phd",
        **opts,
    )
    flopy.mf6.ModflowGwfrch(gwf, stress_period_data={0: [[(0, 0, 0), q]]})
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: [[(nlay - 1, 0, 0), zc[-1]]]})
    flopy.mf6.ModflowGwfoc(
        gwf, head_filerecord=f"{name}.hds", saverecord=[("HEAD", "LAST")]
    )
    return sim, None


def check_output(idx, test):
    name = cases[idx]
    ws = test.workspace
    phd = flopy.utils.HeadFile(ws / f"{name}.phd", text="PRESSHEAD")
    sat = flopy.utils.HeadFile(ws / f"{name}.sat", text="SATURATION")
    psi = phd.get_data().ravel()
    s = sat.get_data().ravel()

    # upper quarter of the column: unit gradient
    upper = slice(0, nlay // 4)
    expected = psi_unit_gradient(idx)
    err = np.abs(psi[upper] - expected).max()
    assert err < 1.0e-4, f"pressure head {psi[upper].mean()} != {expected} ({err})"

    if haevgn[idx] is not None:
        se = se_vg(psi[upper]) / se_vg(haevgn[idx])
        s_expected = satres + (1.0 - satres) * se
        assert np.allclose(s[upper], s_expected, atol=1.0e-10)
        # saturated between the water table and the air-entry head
        near_wt = (psi > haevgn[idx]) & (psi < 0.0)
        assert np.all(s[near_wt] == 1.0)

    # the variants must actually differ from the plain curve
    if idx > 0:
        assert abs(expected - psi_unit_gradient(0)) > 0.1


@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        targets=targets,
    )
    test.run()

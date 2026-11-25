"""
Compare the three VSF soil models (Brooks-Corey, Haverkamp, Van Genuchten)
in one-dimensional columns.

Soils: the Haverkamp and Van Genuchten sands of Celia et al. (1990), and a
Brooks-Corey description of the Van Genuchten sand (air entry -25 cm,
lambda = 0.75). Units are cm and s.

test_infiltration_column
    Infiltration from a wet top boundary into a dry 100 cm column, once per
    soil model, with the Newton formulation. Checks the volume balance, that
    the front has advanced, and that the saturation written by MODFLOW equals
    the soil function evaluated in Python at the simulated pressure head,
    in every cell and at every time step.

test_unit_gradient
    A steady downward flux q = 0.5 Ks through a tall column with the water
    table at the bottom. Far above the water table the gradient is one, so
    kr(psi) = q / Ks there; the simulated pressure head is compared with the
    analytical inverse of each kr function.

With --plot, test_infiltration_column writes soil_functions.png (retention,
relative permeability and moisture capacity of the three models, with the
simulated points) and infiltration_profiles.png (pressure head and water
content profiles of the three models) to the test workspace. Use --keep to
keep the workspace.
"""

import flopy
import numpy as np
import pytest
from gwf_test_utils import get_balance_error, get_vsf_soil_data

hvk = get_vsf_soil_data("Celia1990-eq10-Haverkamp")
vgn = get_vsf_soil_data("Celia1990-eq13-VanGenuchten")

SOILS = {
    "brookscorey": dict(
        porosity=vgn["porosity"],
        satres=vgn["satres"],
        ks=0.00922,
        alpha=0.04,
        lam=0.75,
    ),
    "haverkamp": dict(
        porosity=hvk["porosity"],
        satres=hvk["satres"],
        ks=0.00944,
        alpha=hvk["alpha"],
        n=hvk["n"],
        beta=hvk["beta"],
        k=hvk["k"],
    ),
    "vangenuchten": dict(
        porosity=vgn["porosity"],
        satres=vgn["satres"],
        ks=0.00922,
        alpha=vgn["alpha"],
        n=vgn["n"],
    ),
}
MODELS = list(SOILS)
LABELS = {"brookscorey": "Brooks-Corey", "haverkamp": "Haverkamp"}
LABELS["vangenuchten"] = "Van Genuchten"


# --- soil functions in Python, mirroring the Fortran implementation ---


def se(model, psi):
    """Effective saturation"""
    p = SOILS[model]
    psi = np.asarray(psi, dtype=float)
    h = np.maximum(-psi, 0.0)
    if model == "brookscorey":
        x = p["alpha"] * h
        return np.where(x <= 1.0, 1.0, np.maximum(x, 1.0) ** -p["lam"])
    if model == "haverkamp":
        return 1.0 / (1.0 + (p["alpha"] * h) ** p["n"])
    m = 1.0 - 1.0 / p["n"]
    return (1.0 + (p["alpha"] * h) ** p["n"]) ** -m


def kr(model, psi):
    """Relative permeability"""
    p = SOILS[model]
    psi = np.asarray(psi, dtype=float)
    h = np.maximum(-psi, 0.0)
    if model == "brookscorey":
        return se(model, psi) ** (3.0 + 2.0 / p["lam"])
    if model == "haverkamp":
        return 1.0 / (1.0 + (p["beta"] * h) ** p["k"])
    m = 1.0 - 1.0 / p["n"]
    s = se(model, psi)
    return np.sqrt(s) * (1.0 - (1.0 - s ** (1.0 / m)) ** m) ** 2


def saturation(model, psi):
    sr = SOILS[model]["satres"]
    return sr + (1.0 - sr) * se(model, psi)


def capacity(model, psi, eps=1.0e-6):
    """Moisture capacity d(theta)/d(psi), by central difference"""
    p = SOILS[model]
    ds = saturation(model, psi + eps) - saturation(model, psi - eps)
    return p["porosity"] * ds / (2.0 * eps)


def psi_kr(model, target):
    """Pressure head at which kr equals target (closed form or bisection)"""
    p = SOILS[model]
    if model == "brookscorey":
        return -(target ** (-1.0 / (2.0 + 3.0 * p["lam"]))) / p["alpha"]
    if model == "haverkamp":
        return -((1.0 / target - 1.0) ** (1.0 / p["k"])) / p["beta"]
    lo, hi = -1.0e4, 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if kr(model, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# --- model set-up ---


def add_vsf(gwf, model, name):
    p = SOILS[model]
    pars = dict(porosity=p["porosity"], satres=p["satres"])
    if model == "brookscorey":
        pars.update(soil_model="brooks-corey", alphabc=p["alpha"], lambdabc=p["lam"])
    elif model == "haverkamp":
        pars.update(
            soil_model="haverkamp",
            alphahvk=p["alpha"],
            nhvk=p["n"],
            betahvk=p["beta"],
            khvk=p["k"],
        )
    else:
        pars.update(soil_model="vangenuchten", alphavgn=p["alpha"], nvgn=p["n"])
    flopy.mf6.ModflowGwfvsf(
        gwf,
        ivarsat=1,
        sat_filerecord=f"{name}.sat",
        phead_filerecord=f"{name}.phd",
        **pars,
    )


def column(ws, name, model, exe, nlay, delz, perioddata, steady, strt_psi):
    """A column of nlay cells of delz cm; returns the simulation and centres"""
    top = 0.0
    botm = top - delz * np.arange(1, nlay + 1)
    zc = botm + 0.5 * delz
    sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name=exe)
    flopy.mf6.ModflowTdis(sim, time_units="SECONDS", perioddata=perioddata)
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, newtonoptions="newton")
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        outer_dvclose=1.0e-9,
        outer_maximum=200,
        inner_maximum=300,
        inner_dvclose=1.0e-10,
        rcloserecord=1.0e-10,
        linear_acceleration="BICGSTAB",
    )
    sim.register_ims_package(ims, [gwf.name])
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=1, ncol=1, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=zc + strt_psi)
    ks = SOILS[model]["ks"]
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=ks, k33=ks)
    if steady:
        flopy.mf6.ModflowGwfsto(gwf, iconvert=0, steady_state={0: True})
    else:
        flopy.mf6.ModflowGwfsto(gwf, iconvert=0, ss=0.0, transient={0: True})
    add_vsf(gwf, model, name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )
    return sim, gwf, zc


def run(sim):
    sim.write_simulation(silent=True)
    success, buff = sim.run_simulation(silent=True, report=True)
    assert success, "\n".join(buff[-20:])


# --- infiltration into a dry column ---

INF_NLAY, INF_DZ = 100, 1.0
INF_PSI0, INF_PSITOP = -500.0, -20.0
INF_PERIOD = (1200.0, 40, 1.0)  # 20 minutes in steps of 30 s


def test_infiltration_column(function_tmpdir, targets, plot):
    results = {}
    for model in MODELS:
        name = f"inf-{model[:3]}"
        ws = function_tmpdir / model
        sim, gwf, zc = column(
            ws,
            name,
            model,
            targets["mf6"],
            INF_NLAY,
            INF_DZ,
            [INF_PERIOD],
            False,
            INF_PSI0,
        )
        flopy.mf6.ModflowGwfchd(
            gwf, stress_period_data=[[(0, 0, 0), zc[0] + INF_PSITOP]]
        )
        run(sim)

        phd = flopy.utils.HeadFile(ws / f"{name}.phd", text="PRESSHEAD")
        sat = flopy.utils.HeadFile(ws / f"{name}.sat", text="SATURATION")
        times = phd.get_times()
        psi = np.array([phd.get_data(totim=t).ravel() for t in times])
        s = np.array([sat.get_data(totim=t).ravel() for t in times])
        results[model] = dict(times=np.array(times), zc=zc, psi=psi, sat=s)

        # MODFLOW's saturation is the soil function at the simulated head
        ds = np.abs(s - saturation(model, psi)).max()
        assert ds < 1.0e-10, f"{model}: saturation differs from soil function {ds}"

        # volume balance
        err = get_balance_error(ws / f"{name}.lst")
        assert abs(err) < 1.0e-3, f"{model}: cumulative balance error {err} %"

        # the front has moved in but not reached the bottom; ahead of it
        # the uniform initial head drains slowly under gravity (< 1e-4)
        s0 = saturation(model, INF_PSI0)
        assert s[-1, 5] > s0 + 0.1, f"{model}: no infiltration front"
        assert s[-1, -1] - s0 < 1.0e-3, f"{model}: front reached the bottom"

    if plot:
        plot_soil_functions(function_tmpdir, results)
        plot_profiles(function_tmpdir, results)


# --- steady unit-gradient column ---

UG_NLAY, UG_DZ = 100, 2.0
UG_KR = 0.5


@pytest.mark.parametrize("model", MODELS)
def test_unit_gradient(model, function_tmpdir, targets):
    name = f"ug-{model[:3]}"
    sim, gwf, zc = column(
        function_tmpdir,
        name,
        model,
        targets["mf6"],
        UG_NLAY,
        UG_DZ,
        [(1.0, 1, 1.0)],
        True,
        -10.0,
    )
    q = UG_KR * SOILS[model]["ks"]
    flopy.mf6.ModflowGwfrch(gwf, stress_period_data={0: [[(0, 0, 0), q]]})
    flopy.mf6.ModflowGwfchd(
        gwf, stress_period_data={0: [[(UG_NLAY - 1, 0, 0), zc[-1]]]}
    )
    run(sim)

    psi = flopy.utils.HeadFile(
        function_tmpdir / f"{name}.phd", text="PRESSHEAD"
    ).get_data()
    psi = psi.ravel()
    expected = psi_kr(model, UG_KR)
    upper = psi[: UG_NLAY // 4]
    err = np.abs(upper - expected).max()
    assert err < 1.0e-4, f"{model}: psi {upper.mean()} != {expected} ({err})"


# --- plots (--plot) ---


def plot_soil_functions(ws, results):
    import matplotlib.pyplot as plt

    h = np.logspace(-1, 4, 400)
    psi = -h
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2), layout="constrained")
    for i, model in enumerate(MODELS):
        color = f"C{i}"
        p = SOILS[model]
        theta = p["porosity"] * saturation(model, psi)
        axs[0].plot(h, theta, color=color, label=LABELS[model])
        r = results[model]
        sim_h = -r["psi"][r["psi"] < 0.0]
        sim_t = p["porosity"] * r["sat"][r["psi"] < 0.0]
        axs[0].plot(sim_h, sim_t, ".", color=color, ms=2, alpha=0.3)
        axs[1].plot(h, kr(model, psi), color=color, label=LABELS[model])
        axs[2].plot(h, capacity(model, psi), color=color, label=LABELS[model])
    axs[0].set_ylabel("water content θ (-)")
    axs[0].set_title("retention (dots: simulated)")
    axs[1].set_ylabel("relative permeability kr (-)")
    axs[1].set_yscale("log")
    axs[1].set_ylim(1.0e-10, 2.0)
    axs[1].set_title("relative permeability")
    axs[2].set_ylabel("moisture capacity dθ/dψ (1/cm)")
    axs[2].set_yscale("log")
    axs[2].set_ylim(1.0e-8, 1.0e-1)
    axs[2].set_title("moisture capacity")
    for ax in axs:
        ax.set_xscale("log")
        ax.set_xlabel("suction -ψ (cm)")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
    fname = ws / "soil_functions.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"wrote {fname}")


def plot_profiles(ws, results):
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(10, 6), sharey=True, layout="constrained")
    styles = ["-", "--", "-.", ":"]
    for i, model in enumerate(MODELS):
        r = results[model]
        p = SOILS[model]
        nt = len(r["times"])
        for j, it in enumerate(np.linspace(nt // 8, nt - 1, 4).astype(int)):
            label = f"{LABELS[model]}, {r['times'][it] / 60:.0f} min"
            axs[0].plot(r["psi"][it], r["zc"], styles[j], color=f"C{i}", label=label)
            axs[1].plot(p["porosity"] * r["sat"][it], r["zc"], styles[j], color=f"C{i}")
    axs[0].set_xlabel("pressure head ψ (cm)")
    axs[0].set_ylabel("elevation (cm)")
    axs[1].set_xlabel("water content θ (-)")
    axs[0].legend(fontsize=7, ncol=1)
    for ax in axs:
        ax.grid(True, alpha=0.3)
    fig.suptitle("Infiltration into a dry column")
    fname = ws / "infiltration_profiles.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"wrote {fname}")

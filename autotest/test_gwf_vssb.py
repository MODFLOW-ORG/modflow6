"""
Test the variably saturated surface boundary (VSSB) package.

VSSB applies precipitation to the top face of surface cells of a Richards (VSF)
model.  Infiltration is limited by the soil through a conductance that is
scaled by the relative permeability, water that cannot infiltrate is stored in
a pond of limited depth, groundwater can exfiltrate into the pond, and water
that cannot be stored leaves the package as runoff.

All models are 1D columns of Richards cells with a constant head at the bottom
cell.  Lengths are in cm and times in seconds.  The balance of the package,

    PRECIPITATION = GWF + RUNOFF + RUNOFF-TO-MVR + STORAGE + EVAPORATION,

is verified with the package observations for every time step.

Cases:

* under:   precipitation below the infiltration capacity; no pond and no runoff,
           and the heads equal those of the same flux applied with RCH
* over:    precipitation above the capacity; the pond fills to its maximum depth
           and the excess leaves as runoff
* nopond:  as over, with a maximum pond depth of zero; all excess is runoff
* picard:  as over, without the Newton formulation
* exfil:   water table above the land surface; the pond fills and then runs off,
           and the steady state satisfies Q = Cs * (zs + dmax - h)
* pulse:   a precipitation pulse followed by a dry period; the pond drains
* steady:  a steady-state period matches the long-time transient solution
* chd:     the surface cell is a constant head cell; all precipitation is runoff
* mvr:     runoff is moved to a stream reach
* ts:      precipitation from a time series, with an auxiliary multiplier
* evap-pond: precipitation ponds and evaporates at the potential rate; all
           evaporation comes from the pond and the soil evaporation is zero
* evap-soil: no precipitation and no pond; the potential evaporation is applied
           to the soil, extracted from the aquifer and reduced by the relative
           permeability
"""

import flopy
import numpy as np
import pandas as pd
import pytest
from framework import TestFramework
from gwf_test_utils import get_balance_error, get_vsf_soil_data

cases = [
    "vssb-under",
    "vssb-over",
    "vssb-nopond",
    "vssb-picard",
    "vssb-exfil",
    "vssb-pulse",
    "vssb-steady",
    "vssb-chd",
    "vssb-mvr",
    "vssb-ts",
    "vssb-evap-pond",
    "vssb-evap-soil",
]

nlay = 20
delz = 1.0  # cm
area = 1.0  # cm2
top = nlay * delz
hk = 0.00944  # cm/s
dmax_default = 0.5  # cm
# conductance between the surface and the center of the top cell at krel = 1
cs_sat = hk * area / (0.5 * delz)
mvr_factor = 0.6  # fraction of the runoff that is moved
obs_names = [
    "precipitation",
    "infiltration",
    "runoff",
    "runoff-to-mvr",
    "pond-storage",
    "evaporation",
    "pond-depth",
]

# case settings: precipitation (cm/s), maximum pond depth (cm), head of the
# constant head cell at the bottom (cm), period length (s), number of steps
settings = {
    "vssb-under": dict(precip=1.0e-4, dmax=dmax_default, perlen=600.0, nstp=60),
    "vssb-over": dict(precip=0.05, dmax=dmax_default, perlen=600.0, nstp=60),
    "vssb-nopond": dict(precip=0.05, dmax=0.0, perlen=600.0, nstp=60),
    "vssb-picard": dict(
        precip=0.05, dmax=dmax_default, perlen=600.0, nstp=60, newton=False
    ),
    "vssb-exfil": dict(
        precip=0.0, dmax=0.1, chd=top + 0.3, perlen=5000.0, nstp=250, strt=top + 0.3
    ),
    "vssb-pulse": dict(
        precip=0.05, dmax=dmax_default, perlen=[100.0, 800.0], nstp=[20, 80]
    ),
    "vssb-steady": dict(precip=0.002, dmax=dmax_default, perlen=1.0, nstp=1),
    "vssb-chd": dict(precip=0.05, dmax=dmax_default, perlen=100.0, nstp=10),
    "vssb-mvr": dict(precip=0.05, dmax=dmax_default, perlen=600.0, nstp=60),
    "vssb-ts": dict(precip=0.0, dmax=dmax_default, perlen=600.0, nstp=60),
    "vssb-evap-pond": dict(
        precip=0.05, evap=0.002, dmax=dmax_default, perlen=600.0, nstp=60
    ),
    "vssb-evap-soil": dict(
        precip=0.0,
        evap=0.002,
        dmax=0.0,
        chd=top - 1.0,
        strt=top - 1.0,
        perlen=600.0,
        nstp=60,
    ),
}


def get_settings(name):
    s = dict(
        precip=0.0,
        evap=0.0,
        dmax=dmax_default,
        chd=0.0,
        strt=0.0,
        perlen=600.0,
        nstp=60,
        newton=True,
    )
    s.update(settings[name])
    s["perlen"] = np.atleast_1d(s["perlen"])
    s["nstp"] = np.atleast_1d(s["nstp"])
    return s


def build_column(
    ws,
    name,
    s,
    vssb=True,
    rch=False,
    steady=False,
    top_chd=False,
    mvr=False,
    ts=False,
    extra=None,
):
    """Build a 1D column of Richards cells with a constant head at the bottom."""
    nper = len(s["perlen"])
    ncol = 2 if mvr else 1
    sim = flopy.mf6.MFSimulation(
        sim_name=name, version="mf6", exe_name="mf6", sim_ws=ws
    )
    flopy.mf6.ModflowTdis(
        sim,
        time_units="SECONDS",
        nper=nper,
        perioddata=[(s["perlen"][i], s["nstp"][i], 1.0) for i in range(nper)],
    )
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="SUMMARY",
        complexity="COMPLEX",
        inner_dvclose=1e-9,
        outer_dvclose=1e-9,
        outer_maximum=300,
        inner_maximum=300,
    )

    gwfname = "gwf"
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=gwfname,
        save_flows=True,
        newtonoptions="newton" if s["newton"] else None,
    )
    sim.register_ims_package(ims, [gwf.name])
    botm = [top - (i + 1) * delz for i in range(nlay)]
    flopy.mf6.ModflowGwfdis(
        gwf, nlay=nlay, nrow=1, ncol=ncol, delr=1.0, delc=1.0, top=top, botm=botm
    )
    flopy.mf6.ModflowGwfic(gwf, strt=s["strt"])
    flopy.mf6.ModflowGwfnpf(gwf, icelltype=0, k=hk, k33=hk)
    flopy.mf6.ModflowGwfsto(
        gwf,
        ss=0.0,
        sy=0.0,
        iconvert=0,
        steady_state={0: True} if steady else None,
        transient=None if steady else {0: True},
    )
    soil = get_vsf_soil_data("Celia1990-eq10-Haverkamp")
    flopy.mf6.ModflowGwfvsf(
        gwf,
        iunsat=1,
        storage_scheme="chord-slope",
        kr_averaging="arithmetic",
        soil_model="Haverkamp",
        porosity=soil["porosity"],
        satres=soil["satres"],
        alphahvk=soil["alpha"],
        nhvk=soil["n"],
        betahvk=soil["beta"],
        khvk=soil["k"],
    )

    chd = [[(nlay - 1, 0, 0), s["chd"]]]
    if top_chd:
        chd.append([(0, 0, 0), s["chd"] + 5.0])
    if mvr:
        chd.append([(nlay - 1, 0, 1), s["chd"]])
    flopy.mf6.ModflowGwfchd(gwf, stress_period_data={0: chd}, pname="CHD-1")

    if rch:
        flopy.mf6.ModflowGwfrch(
            gwf, stress_period_data={0: [[(0, 0, 0), s["precip"]]]}, pname="RCH-1"
        )

    if vssb:
        # flopy numbers surface cells from zero
        options = dict(
            pname="VSSB-1",
            print_flows=True,
            save_flows=True,
            pond_filerecord=f"{gwfname}.pond",
            observations={
                f"{gwfname}.vssb.obs.csv": [
                    (n, n, 1) for n in obs_names if mvr or n != "runoff-to-mvr"
                ]
            },
        )
        perioddata = {0: [[0, s["precip"], s["evap"]]]}
        if len(s["perlen"]) > 1:
            perioddata[1] = [[0, 0.0, s["evap"]]]
        if ts:
            options["timeseries"] = {
                "filename": f"{gwfname}.vssb.ts",
                "timeseries": [(0.0, 0.0), (100.0, 0.05), (200.0, 0.0), (1e5, 0.0)],
                "time_series_namerecord": "rain",
                "interpolation_methodrecord": "stepwise",
            }
            options["auxiliary"] = ["mult"]
            options["auxmultname"] = "mult"
            perioddata = {0: [[0, "rain", s["evap"], 2.0]]}
        if mvr:
            options["mover"] = True
        flopy.mf6.ModflowGwfvssb(
            gwf,
            nsurfcells=1,
            packagedata=[[0, (0, 0, 0), s["dmax"]]],
            perioddata=perioddata,
            **options,
        )

    if mvr:
        flopy.mf6.ModflowGwfsfr(
            gwf,
            nreaches=1,
            packagedata=[
                [0, (0, 0, 1), 1.0, 1.0, 1e-3, top - 0.2, 0.5, 0.0, 0.03, 0, 0, 0]
            ],
            connectiondata=[[0]],
            pname="SFR-1",
            mover=True,
            observations={
                f"{gwfname}.sfr.obs.csv": [("frommvr", "from-mvr", 1)],
            },
        )
        flopy.mf6.ModflowGwfmvr(
            gwf,
            maxpackages=2,
            packages=[("VSSB-1",), ("SFR-1",)],
            maxmvr=1,
            perioddata=[["VSSB-1", 0, "SFR-1", 0, "FACTOR", mvr_factor]],
        )

    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{gwfname}.cbc",
        head_filerecord=f"{gwfname}.hds",
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("BUDGET", "LAST")],
    )
    if extra is not None:
        extra(gwf)
    return sim


def read_obs(ws, fname="gwf.vssb.obs.csv"):
    return pd.read_csv(ws / fname).set_index("time")


def check_balance(obs, mvr=False):
    """The package budget closes for every time step.

    This obs-based identity holds when there is no bare-soil evaporation (which
    is extracted from the aquifer and so is folded into the GWF term rather than
    into the pond exchange reported by the infiltration observation).  Cases with
    soil evaporation are verified through the model listing-file balance instead.
    """
    rhs = obs["INFILTRATION"] - obs["RUNOFF"] - obs["POND-STORAGE"] - obs["EVAPORATION"]
    if mvr:
        rhs = rhs - obs["RUNOFF-TO-MVR"]
    assert np.allclose(obs["PRECIPITATION"], rhs, atol=1e-10, rtol=0.0), (
        "package budget not closed"
    )


def max_outer_iterations(ws):
    """Largest number of outer iterations needed by a time step."""
    counts = []
    for line in open(ws / "mfsim.lst"):
        if "CALLS TO NUMERICAL SOLUTION" in line:
            counts.append(int(line.split()[0]))
    return max(counts)


def build_models(idx, test):
    name = cases[idx]
    s = get_settings(name)
    ws = test.workspace
    kwargs = {
        "vssb-steady": dict(steady=True),
        "vssb-chd": dict(top_chd=True),
        "vssb-mvr": dict(mvr=True),
        "vssb-ts": dict(ts=True),
    }.get(name, {})
    sim = build_column(ws, name, s, **kwargs)

    # reference models run in their own workspaces
    if name == "vssb-under":
        build_column(ws / "rch", name, s, vssb=False, rch=True).write_simulation()
    if name == "vssb-steady":
        # transient run to a long time
        s2 = dict(s, perlen=np.array([5.0e5]), nstp=np.array([500]))
        s2["perlen"] = np.array([5.0e5])
        sim2 = build_column(ws / "tr", name, s2)
        sim2.write_simulation()
    return sim, None


def run_reference(ws, name="vssb"):
    sim = flopy.mf6.MFSimulation.load(sim_ws=str(ws), verbosity_level=0)
    sim.exe_name = "mf6"
    success, _ = sim.run_simulation(silent=True)
    assert success, f"reference model in {ws} failed"


def check_output(idx, test):
    name = cases[idx]
    s = get_settings(name)
    ws = test.workspace
    obs = read_obs(ws)
    heads = flopy.utils.HeadFile(ws / "gwf.hds").get_data(idx=-1).flatten()
    mvr = name == "vssb-mvr"
    # the obs-based package balance excludes soil evaporation (folded into GWF);
    # the soil-evaporation case is verified through the listing-file balance
    if name != "vssb-evap-soil":
        check_balance(obs, mvr=mvr)
    qp = s["precip"]

    # the model budget closes and the solution converges in a bounded number
    # of outer iterations
    assert abs(get_balance_error(ws / "gwf.lst")) < 1e-3
    assert max_outer_iterations(ws) < 50

    # the pond depth file agrees with the observation
    pond = flopy.utils.HeadFile(ws / "gwf.pond", text="POND DEPTH")
    assert np.isclose(pond.get_data(idx=-1)[0, 0, 0], obs["POND-DEPTH"].iloc[-1])

    if name == "vssb-under":
        assert (obs["POND-DEPTH"] == 0.0).all(), "no pond expected"
        assert (obs["RUNOFF"] == 0.0).all(), "no runoff expected"
        assert np.allclose(obs["INFILTRATION"], qp * area, rtol=1e-8)
        # the same flux applied with RCH
        run_reference(ws / "rch")
        href = flopy.utils.HeadFile(ws / "rch" / "gwf.hds").get_data(idx=-1).flatten()
        assert np.allclose(heads, href, atol=1e-6), "heads differ from RCH"

    elif name in ("vssb-over", "vssb-picard"):
        assert np.isclose(obs["POND-DEPTH"].iloc[-1], s["dmax"])
        assert obs["RUNOFF"].iloc[-1] < 0.0, "runoff expected"
        if name == "vssb-over":
            # Newton and Picard agree
            ws_picard = ws / "picard"
            build_column(ws_picard, name, dict(s, newton=False)).write_simulation()
            run_reference(ws_picard)
            href = (
                flopy.utils.HeadFile(ws_picard / "gwf.hds").get_data(idx=-1).flatten()
            )
            assert np.allclose(heads, href, atol=1e-6), "Newton and Picard differ"
            oref = read_obs(ws_picard)
            assert np.allclose(obs, oref, rtol=1e-4, atol=1e-8)
        assert np.isclose(
            -obs["RUNOFF"].iloc[-1], qp * area - obs["INFILTRATION"].iloc[-1]
        )
        # the infiltration is limited by the capacity of the soil
        assert 0.0 < obs["INFILTRATION"].iloc[-1] < qp * area

    elif name == "vssb-nopond":
        assert (obs["POND-DEPTH"] == 0.0).all()
        assert np.allclose(
            -obs["RUNOFF"], qp * area - obs["INFILTRATION"], atol=1e-12
        ), "all water that does not infiltrate must be runoff"
        assert (obs["POND-STORAGE"] == 0.0).all()

    elif name == "vssb-exfil":
        # steady state: Q = Cs * (zs + dmax - h), exfiltration is negative
        h_top = heads[0]
        q_expected = cs_sat * (top + s["dmax"] - h_top)
        assert q_expected < 0.0
        assert np.isclose(obs["INFILTRATION"].iloc[-1], q_expected, rtol=1e-4)
        assert np.isclose(obs["POND-DEPTH"].iloc[-1], s["dmax"])
        assert np.isclose(-obs["RUNOFF"].iloc[-1], -q_expected, rtol=1e-4)

    elif name == "vssb-pulse":
        assert obs["POND-DEPTH"].max() > 0.0
        assert obs["POND-DEPTH"].iloc[-1] == 0.0, "pond must drain"
        assert (obs["POND-DEPTH"] <= s["dmax"] + 1e-12).all()
        # no chatter: the pond depth does not increase after the pulse
        after = obs.loc[obs.index > s["perlen"][0], "POND-DEPTH"].to_numpy()
        assert (np.diff(after) <= 1e-12).all()

    elif name == "vssb-steady":
        run_reference(ws / "tr")
        otr = read_obs(ws / "tr")
        assert np.isclose(
            obs["INFILTRATION"].iloc[-1], otr["INFILTRATION"].iloc[-1], rtol=1e-4
        )
        htr = flopy.utils.HeadFile(ws / "tr" / "gwf.hds").get_data(idx=-1).flatten()
        assert np.allclose(heads, htr, atol=1e-3)

    elif name == "vssb-chd":
        assert (obs["INFILTRATION"] == 0.0).all()
        assert (obs["POND-DEPTH"] == 0.0).all()
        assert np.allclose(-obs["RUNOFF"], qp * area)

    elif name == "vssb-mvr":
        sfr = pd.read_csv(ws / "gwf.sfr.obs.csv").set_index("time")
        assert (obs["RUNOFF-TO-MVR"] < 0.0).any()
        assert np.allclose(sfr["FROMMVR"], -obs["RUNOFF-TO-MVR"], rtol=1e-6)
        # the part that is not moved leaves the model as runoff
        total = obs["RUNOFF"] + obs["RUNOFF-TO-MVR"]
        assert (obs["RUNOFF"] < 0.0).any()
        assert np.allclose(obs["RUNOFF-TO-MVR"], mvr_factor * total, rtol=1e-6)

    elif name == "vssb-ts":
        # the multiplier of 2 doubles the rainfall during the pulse
        assert np.isclose(obs["PRECIPITATION"].max(), 2.0 * 0.05 * area)
        assert (obs["PRECIPITATION"].iloc[25:] == 0.0).all()

    elif name == "vssb-evap-pond":
        # the pond fills and open water evaporates at the potential rate; while
        # the pond holds water there is no bare-soil evaporation, so the obs
        # package balance (which excludes soil evaporation) still closes
        epot = s["evap"] * area
        assert np.isclose(obs["POND-DEPTH"].iloc[-1], s["dmax"])
        assert (obs["EVAPORATION"] <= 0.0).all()
        assert np.isclose(-obs["EVAPORATION"].iloc[-1], epot, rtol=1e-6)

    elif name == "vssb-evap-soil":
        # no pond and no precipitation: all evaporation is bare-soil
        # evaporation, extracted from the aquifer and reduced by the relative
        # permeability, so it is positive but no larger than the potential rate
        epot = s["evap"] * area
        e_soil = -obs["EVAPORATION"].iloc[-1]
        assert (obs["POND-DEPTH"] == 0.0).all(), "no pond expected"
        assert (obs["RUNOFF"] == 0.0).all(), "no runoff expected"
        assert 0.0 < e_soil <= epot + 1e-12, "soil evaporation reduced by krel"
        # the aquifer supplies the evaporation, so the top head drops
        assert heads[0] < s["strt"]


def run_expecting_error(ws, sim, message):
    sim.write_simulation(silent=True)
    success, _ = sim.run_simulation(silent=True)
    assert not success, "the simulation should have failed"
    lst = (ws / "mfsim.lst").read_text()
    assert message in " ".join(lst.split()), f"error message not found: {message}"


def vssp_overlap(ihc):
    def add(gwf):
        flopy.mf6.ModflowGwfvssp(
            gwf,
            maxbound=1,
            stress_period_data={0: [[(0, 0, 0), ihc, 0.5, 1.0, 0.0]]},
            pname="VSSP-1",
        )

    return add


# name, keyword arguments of the model, message in the listing file
validation_cases = {
    "no-vsf": (
        dict(vsf="remove"),
        "The VSSB package requires the VSF package",
    ),
    "duplicate-cell": (
        dict(packagedata=[[0, (0, 0, 0), 0.1], [1, (0, 0, 0), 0.1]], nsurfcells=2),
        "is specified more than once in the VSSB PACKAGEDATA block",
    ),
    "active-above": (
        dict(cellid=(1, 0, 0)),
        "has an active cell directly above it",
    ),
    "vssp-overlap": (
        dict(extra=vssp_overlap(0)),
        "is a VSSB surface cell and cannot be listed with IHC = 0",
    ),
    "negative-precipitation": (
        dict(precip=-0.01),
        "PRECIPITATION for surface cell 1 must be greater than or equal to 0",
    ),
    "negative-evaporation": (
        dict(evap=-0.01),
        "EVAPORATION for surface cell 1 must be greater than or equal to 0",
    ),
    "negative-pond-depth": (
        dict(dmax=-0.1),
        "MAX_POND_DEPTH for surface cell 1 must be greater than or equal to 0",
    ),
}


@pytest.mark.parametrize("name", validation_cases)
def test_validation(name, function_tmpdir):
    """Invalid input is reported at read time."""
    kwargs, message = validation_cases[name]
    s = get_settings("vssb-over")
    ws = function_tmpdir
    sim = build_column(ws, "vssb-" + name, s, extra=kwargs.get("extra"))
    gwf = sim.get_model("gwf")
    if "vsf" in kwargs:
        gwf.remove_package("vsf")
    vssb = gwf.get_package("vssb-1")
    if "packagedata" in kwargs:
        vssb.nsurfcells = kwargs["nsurfcells"]
        vssb.packagedata.set_data(kwargs["packagedata"])
        vssb.perioddata.set_data(
            {0: [[0, s["precip"], s["evap"]], [1, s["precip"], s["evap"]]]}
        )
    if "cellid" in kwargs:
        vssb.packagedata.set_data([[0, kwargs["cellid"], s["dmax"]]])
    if "precip" in kwargs:
        vssb.perioddata.set_data({0: [[0, kwargs["precip"], s["evap"]]]})
    if "evap" in kwargs:
        vssb.perioddata.set_data({0: [[0, s["precip"], kwargs["evap"]]]})
    if "dmax" in kwargs:
        vssb.packagedata.set_data([[0, (0, 0, 0), kwargs["dmax"]]])
    run_expecting_error(ws, sim, message)


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

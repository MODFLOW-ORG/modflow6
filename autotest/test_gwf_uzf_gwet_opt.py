"""
A test of the groundwater ET options in UZF. An issue was reported (#1944)
that GWET was showing up as in inflow in the budget table.  This test employs
a DISU grid.

The test uses 2 GWF models in a single mf6 simulation. Both models use
a simple 2 x 2 x 2 DISU grid.  One constant head cell in the lower layer
is meant to keep the water table in layer 1. One GWF model use the classic
LINEAR_GWET approximation for calculating GWET after accounting for UZET.
The 2nd model uses SQUARE_GWET option.  The SQUARE_GWET option is not
head-dependent; it applies the pET rate over the entire extinction depth
except for a small smoothing interval at the bottom of the extinction depth.
The upper layer is 2.5 m thick and the extinction depth is 2.25.  The
constant head boundary seeks to keep the water table at roughly 1 m
depth.

         Uses LINEAR_GWET option           Uses SQUARE_GWET option
         -----------------------           -----------------------

             +-------+-------+                  +-------+-------+
            /       /       /|                 /       /       /|
           /   1   /   2   / |                /   1   /   2   / |
          /       /       /  |               /       /       /  |
         +-------+-------+ 2 +              +-------+-------+ 2 +
        /       /       /|  /|             /       /       /|  /|
       /   3   /   4   / | / |            /   3   /   4   / | / |
      /       /       /  |/  |           /       /       /  |/  |
     +-------+-------+ 4 + 6 +          +-------+-------+ 4 + 6 +
     |       |       |  /|  /           |       |       |  /|  /
     |   3   |   4   | / | /            |   3   |   4   | / | /
     |       |       |/  |/             |       |       |/  |/
     +-------+-------+ 8 +              +-------+-------+ 8 +
     |       |       |  /               |       |       |  /
     |   7   |   8   | /                |   7   |   8   | /
     |       |       |/                 |       |       |/
     +-------+-------+                  +-------+-------+

"""

# Imports

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

cases = ["uzf-gwet"]

# Model units
length_units = "meters"
time_units = "days"

nper = 2
tdis_rc = []
for _ in range(nper):
    tdis_rc.append((1.0, 1, 1.0))

# spatial discretization data
nlay, nrow, ncol = 2, 2, 2
delr, delc = 1.0, 1.0
top = 5.0
topm = [5.0, 2.5]
topnp = np.array([[topm[0]] * 4 + [topm[1]] * 4]).squeeze()
botm = [2.5, 0.0]
botnp = np.array([[botm[0]] * 4 + [botm[1]] * 4]).squeeze()
strt = 4.0
idomain = np.ones((nlay, nrow, ncol))

# set hk
hk = 100.0  # m/d
sy = 0.2

# solver options
hclose, rclose = 1e-9, 1e-3


def flatten(xss):
    return [x for xs in xss for x in xs]


def filter_nodes(xvc, yvc, iv, xv, yv):
    # Check if iv is empty
    if iv:
        vert_id = max(iv)

        # Add nodes represented by xvc and yvc if not already contained in xv, yv
        for new_pair in zip(xvc, yvc):
            found = False
            for old_pair in zip(xv, yv):
                if old_pair[0] == new_pair[0] and old_pair[1] == new_pair[1]:
                    found = True
                    break

            if not found:
                # if not already present, add the vertex
                iv.append(vert_id + 1)
                vert_id += 1
                xv.append(new_pair[0])
                yv.append(new_pair[1])

    else:
        iv.extend(list(range(0, 4)))
        xv.extend(xvc)
        yv.extend(yvc)

    return iv, xv, yv


def populate_linked_IDs(xvc, yvc, xv, yv):
    # Ensure 4 items contained in the passed list of nodes to get positional IDs for
    assert len(xvc) == 4, "Number of processed vertices is off"

    vertices = []
    for pair in zip(xvc, yvc):  # xvc should have 4 items
        for i, existing_pair in enumerate(zip(xv, yv)):
            if existing_pair[0] == pair[0] and existing_pair[1] == pair[1]:
                vertices.append(i)
                break

    assert len(vertices) == 4, "Number of vertices should be 4"
    return vertices


def buildout_vertex_locations():
    n = -1
    cell_vert_lkup = {}
    # Only need 1 layer's worth of vertices since they can be repeated
    iv = []
    xv = []
    yv = []

    ly = 0
    for i in np.arange(nrow):
        for j in np.arange(ncol):
            # There are two X locations (leverages the fact that delr = 1.0)
            vert_left = float(j) * delr
            vert_right = vert_left + delr
            # There are two Y locations (only 1 row with unit width)
            vert_back = (nrow - i) * delc
            vert_front = vert_back - delc

            # Define each vertex for the cell in the current column
            if ly == 0:
                # left, back
                xvc = [vert_left]
                yvc = [vert_back]
                # right, back
                xvc.append(vert_right)
                yvc.append(vert_back)
                # right, front
                xvc.append(vert_right)
                yvc.append(vert_front)
                # left, front
                xvc.append(vert_left)
                yvc.append(vert_front)

            # Second, only keep vertices that don't already appear in the
            # respective lists
            iv, xv, yv = filter_nodes(xvc, yvc, iv, xv, yv)

            # Store dictionary entry linking cell ID with its vertices
            vertices = populate_linked_IDs(xvc, yvc, xv, yv)
            n += 1
            cell_vert_lkup[n] = vertices

    # Now loop for the remaining layers
    for i in np.arange(nrow):
        for j in np.arange(ncol):
            cell_id = (i * ncol) + j
            # Only need to find the top layer's vertices once
            verts_to_use = cell_vert_lkup[cell_id]
            for ly in np.arange(1, nlay):
                n = (ly * nrow * ncol) + (i * ncol) + j
                cell_vert_lkup[n] = verts_to_use

    return iv, xv, yv, cell_vert_lkup


def set_connectiondata(n, lay, top, back, left, right, front, bottom):
    # Instantiate empty lists
    jas = [n]
    ihc = [lay]
    cl12 = [n]
    hwva = [n]
    angldeg = [360.0]

    # Calculate half-cell thickness up front for vertical connections
    if lay == 0:
        cl12_val = (topm[lay] - botm[lay]) / 2
    else:
        cl12_val = (botm[lay - 1] - botm[lay]) / 2

    if top:
        jas.append(n - (ncol * nrow))
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, area is 1.0m x 1.0m
        angldeg.append(0.0)  # placeholder only for vertical connections

    if back:
        jas.append(n - ncol)
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a column
        hwva.append(delr)  # the width perpendicular to the connection
        # for horizontal connection, 90.0 deg points in the positive y direction
        angldeg.append(90.0)

    if left:
        jas.append(n - 1)  # left
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delr / 2)  # half the cell width along a horizontal connection (row)
        # for horizontal connection, value of hwva is width along a row
        hwva.append(delr)
        # for horizontal connection, value is 180.0 along negative x-axis
        angldeg.append(180.0)

    if right:
        jas.append(n + 1)  # right
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delr / 2)  # half the cell width along a horizontal connection row
        # for horizontal connection, value of hwva is width along a row
        hwva.append(delc)
        # for horizontal connection, value is 0.0 along positive x-axis
        angldeg.append(0.0)

    if front:
        jas.append(n + ncol)
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a column
        hwva.append(delr)  # the width perpendicular to the connection
        # for horizontal connection, 90.0 deg points in the positive y direction
        angldeg.append(2700.0)

    if bottom:
        jas.append(n + (nrow * ncol))  # below
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, value of hwva is area
        angldeg.append(0.0)  # placeholder only for vertical connections

    return jas, ihc, cl12, hwva, angldeg


def get_conndat(n, lay, row, col):
    top = False
    back = False
    left = False
    right = False
    front = False
    bottom = False

    # iac is the number of connections (plus 1) for each cell
    iac = 1  # start with the "plus 1"
    iac += 3  # all cells have exactly 3 neighbors for a 2 x 2 x 2 grid

    if lay == 0:
        bottom = True
    elif lay == 1:
        top = True

    if row == 0:
        front = True
    elif row == 1:
        back = True

    if col == 0:
        right = True
    elif col == 1:
        left = True

    (jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals) = set_connectiondata(
        n, lay, top, back, left, right, front, bottom
    )

    return iac, jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals


def append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup):
    # Get the vertex IDs for the current cell

    # Start with calculating the x location of the cell
    col_id = n % ncol
    cell_x = delr / 2 + col_id * delr

    # The y location of the cell is fixed at delc/2
    add_inc = 0.0
    if n in [0, 1, 4, 5]:
        add_inc = 1.0
    cell_y = delc / 2 + add_inc

    # Get associated vertices
    vertices = cell_vert_lkup[n]

    # Every cell will have 4 vertices
    new_cell2d = [[n, cell_x, cell_y, 4], vertices]

    # Flatten
    new_cell2d = flatten(new_cell2d)

    return new_cell2d


iv_lst = []
xv_lst = []
yv_lst = []
cell_2d_lst = []

# Build out VERTICES block (only needs to happen once)
iv, xv_lst, yv_lst, cell_vert_lkup = buildout_vertex_locations()

vertices = []
for i in np.arange(len(iv)):
    vertices.append([iv[i], xv_lst[i], yv_lst[i]])


iac_lst = []
ja_lst = []
ihc_lst = []
cl12_lst = []
hwva_lst = []
angldeg_lst = []

# Cycle through each layer, row, and column.
for lay in np.arange(nlay):
    for row in np.arange(nrow):
        for col in np.arange(ncol):
            n = (lay * nrow * ncol) + (row * ncol) + col  # n will be zero based

            # Values for CONNECTIONDATA block
            (iac, ja_cell, ihc_cell, cl12_cell, hwva_cell, angldeg_cell) = get_conndat(
                n, lay, row, col
            )

            # accumulate connection information in lists
            iac_lst.append(iac)
            ja_lst.append(ja_cell)
            ihc_lst.append(ihc_cell)
            cl12_lst.append(cl12_cell)
            hwva_lst.append(hwva_cell)
            angldeg_lst.append(angldeg_cell)

            # Add Cell2D information
            cell_2d = append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup)
            cell_2d_lst.append(cell_2d)


iacnp = np.array(iac_lst)
janp = np.array(flatten(ja_lst))
ihcnp = np.array(flatten(ihc_lst))
cl12np = np.array(flatten(cl12_lst))
hwvanp = np.array(flatten(hwva_lst))
angldegnp = np.array(flatten(angldeg_lst))

# Establish flow boundary conditions.
chd_spd_block = []
chd_spd_block.append([4, 4.0])
chd_spd = {0: chd_spd_block}


def build_gwf_model(sim, name, linear_et_flg=True):
    # Instantiate flow model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        newtonoptions="NEWTON",
    )

    # Instantiate solver
    imsgwf = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
        linear_acceleration="BICGSTAB",
        filename=f"{name}.ims",
    )
    sim.register_ims_package(imsgwf, [gwf.name])

    # Instantiate unstructured grid
    flopy.mf6.ModflowGwfdisu(
        gwf,
        length_units=length_units,
        nogrb=True,
        nodes=nlay * nrow * ncol,
        nja=len(janp),
        nvert=len(iv),
        top=topnp,
        bot=botnp,
        area=delr * delc,
        idomain=idomain.flatten(),
        iac=iacnp,
        ja=janp,
        ihc=ihcnp,
        cl12=cl12np,
        hwva=hwvanp,
        angldegx=angldegnp,
        vertices=vertices,
        cell2d=cell_2d_lst,
        pname="DISU",
    )

    # Instantiate the node property flow package
    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=1,  # >0 means saturated thickness varies with computed head
        k=hk,
        k33=hk,
    )

    # Instantiate the storage package
    flopy.mf6.ModflowGwfsto(
        gwf,
        sy=sy,
        ss=1e-5,
        transient={0: True},
    )

    # Instantiate the initial conditions
    flopy.mf6.ModflowGwfic(gwf, strt=strt)

    # Instantiate constant heads for keeping the heads from rising
    flopy.mf6.ModflowGwfchd(
        gwf,
        stress_period_data=chd_spd,
        pname="CHD",
    )

    # UZF package
    # <ifno> <cellid> <landflag> <ivertcon> <surfdep> <vks> <thtr> <thts> <thti> <eps>
    uzf_pkdat = [
        [0, 0, 1, -1, 1e-3, 1, 0.05, 0.25, 0.05, 4, "uzf01"],
        [1, 1, 1, -1, 1e-3, 1, 0.05, 0.25, 0.05, 4, "uzf02"],
        [2, 2, 1, -1, 1e-3, 1, 0.05, 0.25, 0.15, 4, "uzf03"],
        [3, 3, 1, -1, 1e-3, 1, 0.05, 0.25, 0.15, 4, "uzf04"],
    ]

    extdp = 2.25
    pet = 0.001
    zero = 0.0
    # <ifno> <finf> <pet> <extdp> <extwc> <ha> <hroot> <rootact>
    uzf_spd = {
        0: [
            [0, 0.00, pet, extdp, 5.0e-02, zero, zero, zero],
            [1, 0.01, pet, extdp, 5.0e-02, zero, zero, zero],
            [2, 0.01, pet, extdp, 5.0e-02, zero, zero, zero],
            [3, 0.01, pet, extdp, 5.0e-02, zero, zero, zero],
        ]
    }

    if linear_et_flg:
        linear_gwet = True
        square_gwet = False
    else:
        linear_gwet = False
        square_gwet = True

    flopy.mf6.ModflowGwfuzf(
        gwf,
        print_flows=True,
        save_flows=True,
        simulate_et=True,
        simulate_gwseep=False,
        linear_gwet=linear_gwet,
        square_gwet=square_gwet,
        unsat_etwc=True,
        boundnames=True,
        ntrailwaves=15,
        nwavesets=40,
        nuzfcells=len(uzf_pkdat),
        packagedata=uzf_pkdat,
        perioddata=uzf_spd,
        budget_filerecord=f"{name}.uzf.bud",
        pname="UZF",
        filename=f"{name}.uzf",
    )

    # Output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        budget_filerecord=f"{name}.cbc",
        head_filerecord=f"{name}.hds",
        headprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    return gwf


def build_models(idx, test):
    name = cases[idx]
    ws = test.workspace

    # Instantiate a MODFLOW 6 simulation
    sim = flopy.mf6.MFSimulation(
        sim_name=name,
        sim_ws=ws,
        exe_name="mf6",
        version="mf6",
    )

    # Instantiate time discretization package
    flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=tdis_rc, time_units=time_units)

    build_gwf_model(sim, name + "-lin", linear_et_flg=True)
    build_gwf_model(sim, name + "-sqr", linear_et_flg=False)

    return sim, None


def check_output(idx, test):
    msg0 = "GWET should always be negative"
    msg1 = "Results from SQUARE_GWET formulation should be larger"
    msg2 = "no water content in the first uzf object, so should be no uzet"
    msg3 = (
        "Given the shallow depth to the water table, the GWET should "
        "be greater than the UZET"
    )
    msg4 = (
        "The square GWET formulation withdrawals more water and should "
        "result in slightly less heads"
    )

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwfname1 = cases[idx] + "-lin"
    gwfname2 = cases[idx] + "-sqr"
    gwf1 = sim.get_model(gwfname1)
    gwf2 = sim.get_model(gwfname2)
    uzf1 = gwf1.get_package("UZF")
    uzf2 = gwf2.get_package("UZF")

    heads_lin = gwf1.output.head().get_alldata()
    heads_sqr = gwf2.output.head().get_alldata()
    gwet1 = gwf1.output.budget().get_data(text="UZF-GWET")
    gwet2 = gwf2.output.budget().get_data(text="UZF-GWET")
    uzet1 = uzf1.output.budget().get_data(text="UZET")
    uzet2 = uzf2.output.budget().get_data(text="UZET")

    # Because the square gwet formulation removes more water from the model
    # than the linear gwet formulation, the heads should be slightly less
    # in the square_gwet formulation.  However, head differences will be
    # relatively small because of the placement of constant head cell in
    # the lower layer and specified K values.
    for tm in np.arange(nper):
        for ifno in np.arange(heads_sqr.shape[3] - 1):
            assert heads_sqr[tm, 0, 0, ifno] <= heads_lin[tm, 0, 0, ifno], msg4

    gwet1np = []
    gwet2np = []
    uzet1np = []
    uzet2np = []
    for tm in np.arange(nper):
        for ifno in np.arange(len(gwet1[tm])):
            gwet1np.append(gwet1[tm][ifno][-1])
            gwet2np.append(gwet2[tm][ifno][-1])
            uzet1np.append(uzet1[tm][ifno][-1])
            uzet2np.append(uzet2[tm][ifno][-1])

    gwet1np = np.array(gwet1np).reshape((nper, len(gwet1[0])))
    gwet2np = np.array(gwet2np).reshape((nper, len(gwet2[0])))
    uzet1np = np.array(uzet1np).reshape((nper, len(uzet1[0])))
    uzet2np = np.array(uzet2np).reshape((nper, len(uzet2[0])))

    for tm in np.arange(nper):
        for ifno in np.arange(gwet1np.shape[1]):
            assert gwet1np[tm][ifno] < 0, msg0
            assert gwet2np[tm][ifno] < 0, msg0
            assert abs(gwet2np[tm, ifno]) > abs(gwet1np[tm, ifno]), msg1

    # Do a few spot checks of UZ
    assert uzet1np[0, 0] == 0, msg2
    assert uzet2np[0, 0] == 0, msg2
    assert uzet1np[1, 0] == 0, msg2
    assert uzet2np[1, 0] == 0, msg2

    # Confirm UZET vals are all 0.0 or negative
    for tm in np.arange(nper):
        for ifno in np.arange(gwet1np.shape[1]):
            assert uzet1np[tm, ifno] <= 0.0
            assert uzet2np[tm, ifno] <= 0.0

    # Lower soil moisture in the first "row" should assure that GWET > UZET
    for tm in np.arange(nper):
        for ifno in [0, 1]:
            assert abs(uzet1np[tm, ifno]) < abs(gwet1np[tm, ifno]), msg3
            assert abs(uzet2np[tm, ifno]) < abs(gwet2np[tm, ifno]), msg3

    # High soil moisture in the second "row" should assure that UZET > GWET
    # for the linear case, but not in the SQUARE_GWET case (remember that
    # pET is first satisfied from UZ soil moisture and the residual pET
    # that is passed to the GWET calculation will be lower)
    for tm in np.arange(nper):
        for ifno in [2, 3]:
            assert abs(uzet1np[tm, ifno]) > abs(gwet1np[tm, ifno]), msg3
            assert abs(uzet2np[tm, ifno]) < abs(gwet2np[tm, ifno]), msg3


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

"""
A test of the keyword NONE for <cellid> in the package data block.
In essence, this tests inactivation of cells which if active, would be connected
to SFR reaches during a simulation. Additionally, this uses a DISU grid.  The
test was developed based on a user having memory access violations for a real-
world model that included DISU and SFR reaches where cells did not exist.

The test has 6 sfr reaches that are not connected to the
groundwater model.


              --->  Direction of flow  --->
  SFR                                                    |<--- 3 cells wide ---->|
   |     +-------+-------+-------+-------+-------+-------+                       +-------+-------+-------+-------+-------+-------+
   v    /       /       /       /       /       /       /|                      /       /       /       /       /       /       /|
 ===================================================================================================================================
      /       /       /       /       /       /       /  |                    /       /       /       /       /       /       /  |
     +-------+-------+-------+-------+-------+-------+   +                   +-------+-------+-------+-------+-------+-------+   +
     |       |       |       |       |       |       |  /|  For these cells, |       |       |       |       |       |       |  /|
     |       |       |       |       |       |       | / |  idomain set to 0 |       |       |       |       |       |       | / |
     |       |       |       |       |       |       |/  |                   |       |       |       |       |       |       |/  |
     +-------+-------+-------+-------+-------+-------+   +                   +-------+-------+-------+-------+-------+-------+   +
     |       |       |       |       |       |       |  /                    |       |       |       |       |       |       |  /
     |       |       |       |       |       |       | /                     |       |       |       |       |       |       | /
     |       |       |       |       |       |       |/                      |       |       |       |       |       |       |/
     +-------+-------+-------+-------+-------+-------+                       +-------+-------+-------+-------+-------+-------+

"""  # noqa

# Imports

import flopy
import numpy as np
import pytest
from framework import TestFramework

# Base simulation and model name and workspace

scheme = "UPSTREAM"
# scheme = "TVD"

cases = ["sfr-inact"]

# Model units
length_units = "meters"
time_units = "days"

nper = 2
tdis_rc = []
for _ in range(nper):
    tdis_rc.append((1.0, 1, 1.0))

# spatial discretization data
nlay, nrow, ncol = 2, 1, 15
delr, delc = 1.0, 1.0
top = 5.0
topm = [5.0, 2.5]
topnp = np.array([[topm[0]] * 15 + [topm[1]] * 15]).squeeze()
botm = [2.5, 0.0]
botnp = np.array([[botm[0]] * 15 + [botm[1]] * 15]).squeeze()
strt = 1.0
idomain = np.ones((nlay, nrow, ncol))
idomain[:, 0, 6:9] = 0

# set hk
hk = 1.0  # m/d
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
    for j in np.arange(ncol):
        # There are two X locations (leverages the fact that delr = 1.0)
        vert_left = float(j)
        vert_right = float(j) + delr
        # There are two Y locations (only 1 row with unit width)
        vert_front = 0.0
        vert_back = 1.0

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

        # Second, only keep vertices that don't already appear in the respective lists
        iv, xv, yv = filter_nodes(xvc, yvc, iv, xv, yv)

        # Store dictionary entry linking cell ID with its vertices
        vertices = populate_linked_IDs(xvc, yvc, xv, yv)
        n += 1
        cell_vert_lkup[n] = vertices

    # Now loop for the remaining 199 layers
    for j in np.arange(ncol):
        # Only need to find the top layer's vertices once
        verts_to_use = cell_vert_lkup[j]
        for ly in np.arange(1, nlay):
            n = (ly * ncol) + j
            cell_vert_lkup[n] = verts_to_use

    return iv, xv, yv, cell_vert_lkup


def set_connectiondata(n, lay, top, left, right, bottom):
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
        jas.append(n - ncol)
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, area is 1.0m x 1.0m
        angldeg.append(0.0)  # placeholder only for vertical connections

    if left:
        jas.append(n - 1)  # left
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a horizontal connection
        # for horizontal connection, value of hwva is width along a row
        hwva.append(delr)
        # for horizontal connection, value is 180.0 along negative x-axis
        angldeg.append(180.0)

    if right:
        jas.append(n + 1)  # right
        ihc.append(1)  # ihc = 1 for horizontal connection
        cl12.append(delc / 2)  # half the cell width along a horizontal connection
        # for horizontal connection, value of hwva is width along a row
        hwva.append(delc)
        # for horizontal connection, value is 0.0 along positive x-axis
        angldeg.append(0.0)

    if bottom:
        jas.append(n + ncol)  # below
        ihc.append(0)  # ihc = 0 for vertical connection
        cl12.append(cl12_val)  # half the cell thickness or a vertical connection
        hwva.append(delr * delc)  # for vertical connection, value of hwva is area
        angldeg.append(0.0)  # placeholder only for vertical connections

    return jas, ihc, cl12, hwva, angldeg


def get_conndat(n, lay, col):
    top = False
    left = False
    right = False
    bottom = False
    num_neigh = 0

    # iac is the number of connections (plus 1) for each cell
    iac = 1  # start with the "plus 1"

    if col < 6 or col > 8:
        # Cells are connected to only one lateral neighbor
        if col not in [0, 5, 9, 14]:
            iac += 2
            num_neigh = 2
            left = True
            right = True
        else:
            iac += 1
            num_neigh = 1
            if col in [0, 9]:
                right = True
            elif col in [5, 14]:
                left = True
    elif col >= 6 and col <= 8:
        # For the inactive cells in the middle of the model
        # there are no horizontal connections
        pass

    # All cells connected to either cell above or below
    if lay == 0:
        bottom = True
    elif lay == 1:
        top = True

    iac += 1
    num_neigh += 1

    (jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals) = set_connectiondata(
        n, lay, top, left, right, bottom
    )

    # If bottom most layer, need to .pop() the last values out of the respective lists
    # This should be done because the bottom connections will always be represented by
    # the final values in the list (at this point in the development anyway)
    # if lay == nlay - 1:
    #    iac -= 1
    #    jas_vals.pop(-1)
    #    ihc_vals.pop(-1)
    #    cl12_vals.pop(-1)
    #    hwva_vals.pop(-1)
    #    angldeg_vals.pop(-1)

    return iac, jas_vals, ihc_vals, cl12_vals, hwva_vals, angldeg_vals


def append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup):
    # Get the vertex IDs for the current cell

    # Start with calculating the x location of the cell
    col_id = n % ncol
    cell_x = delr / 2 + col_id * delr
    # The y location of the cell is fixed at delc/2
    cell_y = delc / 2

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

# Cycle through each layer and column.
for lay in np.arange(nlay):
    for col in np.arange(ncol):
        n = lay * ncol + col  # n will be zero based

        # Values for CONNECTIONDATA block
        (iac, ja_cell, ihc_cell, cl12_cell, hwva_cell, angldeg_cell) = get_conndat(
            n, lay, col
        )

        # accumulate connection information in lists
        iac_lst.append(iac)
        ja_lst.append(ja_cell)
        ihc_lst.append(ihc_cell)
        cl12_lst.append(cl12_cell)
        hwva_lst.append(hwva_cell)
        angldeg_lst.append(angldeg_cell)

        # Add Cell2D information
        # An example of what cell2d information should look like.
        # [997, 0.09539255, 0.022023124999999998, 4, 997, 1177, 1176, 996],
        # [998, 0.09456585, 0.025338825, 4, 998, 1178, 1177, 997],
        # [999, 0.093623925, 0.028623675, 4, 999, 1179, 1178, 998],
        cell_2d = append_cell2d(n, xv_lst, yv_lst, cell_vert_lkup)
        cell_2d_lst.append(cell_2d)


iacnp = np.array(iac_lst)
janp = np.array(flatten(ja_lst))
ihcnp = np.array(flatten(ihc_lst))
cl12np = np.array(flatten(cl12_lst))
hwvanp = np.array(flatten(hwva_lst))
angldegnp = np.array(flatten(angldeg_lst))

# Establish flow boundary conditions.
chd_spd_left_block = []
chd_spd_right_block = []

chd_spd_left_block.append([20, 1.0, 2.0])
chd_spd_right_block.append([29, 1.0, 2.0])

chd_spd_left = {0: chd_spd_left_block}
chd_spd_right = {0: chd_spd_right_block}


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
        nodes=nlay * ncol * nrow,
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
        auxiliary="CONCENTRATION",
        stress_period_data=chd_spd_left,
        pname="WEL-left",
        filename=f"{name}.wel-left",
    )
    flopy.mf6.ModflowGwfchd(
        gwf,
        auxiliary="CONCENTRATION",
        stress_period_data=chd_spd_right,
        pname="WEL-right",
        filename=f"{name}.wel-right",
    )

    # SFR data
    nreaches = 15

    # <ifno> <cellid(ncelldim)> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> ...
    #        <man> <ncon> <ustrf> <ndv>
    rhk = 0.05
    package_data = [
        (0, 0, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 1, 1.0, 0),
        (1, 1, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (2, 2, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (3, 3, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (4, 4, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (5, 5, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (6, -1, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (7, -1, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (8, -1, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (9, 9, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (10, 10, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (11, 11, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (12, 12, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (13, 13, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 2, 1.0, 0),
        (14, 14, delr, 0.25, 1e-3, top - 0.25, 0.25, rhk, 0.02, 1, 1.0, 0),
    ]

    connection_data = [
        (0, -1),
        (1, 0, -2),
        (2, 1, -3),
        (3, 2, -4),
        (4, 3, -5),
        (5, 4, -6),
        (6, 5, -7),
        (7, 6, -8),
        (8, 7, -9),
        (9, 8, -10),
        (10, 9, -11),
        (11, 10, -12),
        (12, 11, -13),
        (13, 12, -14),
        (14, 13),
    ]

    sfr_spd = {
        0: [
            (0, "inflow", 0.1),
        ],
    }

    sfr_obs = {
        f"{name}.sfr.obs.csv": [
            ("gwf", "sfr", (0,)),
            ("outflow", "ext-inflow", (0,)),
            ("depth", "depth", (0,)),
            ("gwf", "sfr", (6,)),
            ("depth", "depth", (6,)),
            ("gwf", "sfr", (7,)),
            ("depth", "depth", (7,)),
            ("gwf", "sfr", (14,)),
            ("outflow", "ext-outflow", (14,)),
            ("depth", "depth", (14,)),
        ],
        "filename": name + ".sfr.obs",
    }

    flopy.mf6.ModflowGwfsfr(
        gwf,
        save_flows=True,
        print_stage=True,
        print_flows=True,
        print_input=True,
        stage_filerecord=f"{name}.sfr.hds",
        budget_filerecord=f"{name}.sfr.cbc",
        length_conversion=1.0,
        time_conversion=86400,
        mover=False,
        nreaches=nreaches,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=sfr_spd,
        observations=sfr_obs,
        pname="SFR-1",
    )

    # Instantiate output control
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=f"{name}.hds",
        budget_filerecord=f"{name}.cbc",
        saverecord=[("head", "all"), ("budget", "all")],
    )

    # ---------------------------
    # Setup transport model
    # ---------------------------

    gwtname = "gwt-" + name
    gwt = flopy.mf6.MFModel(
        sim,
        model_type="gwt6",
        modelname=gwtname,
        model_nam_file=f"{gwtname}.nam",
    )

    # Instantiate solver
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        outer_dvclose=1e-5,
        inner_dvclose=1e-6,
        linear_acceleration="BICGSTAB",
        filename=f"{gwtname}.ims",
    )
    sim.register_ims_package(imsgwt, [gwt.name])

    # Instantiate same unstructured grid as was used in GWF
    flopy.mf6.ModflowGwtdisu(
        gwt,
        nogrb=True,
        nodes=nlay * ncol * nrow,
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

    # Instantiate initial concentrations
    flopy.mf6.ModflowGwtic(
        gwt,
        strt=1.0,
        filename=f"{gwtname}.ic",
    )

    # Instantiate advection package
    flopy.mf6.ModflowGwtadv(gwt, scheme=scheme, filename=f"{gwtname}.adv")

    # Instantiate mass storage and transfer package
    porosity = sy
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity, filename=f"{gwtname}.sto")

    # Instantiate source/sink mixing package
    sourcerecarray = [
        ("WEL-left", "AUX", "CONCENTRATION"),
        ("WEL-right", "AUX", "CONCENTRATION"),
    ]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray, filename=f"{gwtname}.ssm")

    # Instantiate streamflow transport package
    sftpackagedata = []
    for irno in range(ncol):
        t = (irno, 9.99)
        sftpackagedata.append(t)

    sftperioddata = [(0, "STATUS", "CONSTANT"), (0, "CONCENTRATION", 9.99)]

    flopy.mf6.modflow.ModflowGwtsft(
        gwt,
        flow_package_name="SFR-1",
        boundnames=False,
        save_flows=True,
        print_input=True,
        print_flows=True,
        print_concentration=True,
        concentration_filerecord=gwtname + ".sft.bin",
        budget_filerecord=gwtname + ".sft.bud",
        packagedata=sftpackagedata,
        reachperioddata=sftperioddata,
        pname="SFT",
    )

    # output control
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord=f"{gwtname}.cbc",
        concentration_filerecord=f"{gwtname}.ucn",
        concentrationprintrecord=[("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")],
        saverecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
        printrecord=[("CONCENTRATION", "ALL"), ("BUDGET", "ALL")],
    )

    # Instantiate GWF-GWT exchange
    flopy.mf6.ModflowGwfgwt(
        sim,
        exgtype="GWF6-GWT6",
        exgmnamea=name,
        exgmnameb=gwtname,
        filename=f"{name}.gwfgwt",
    )

    return sim, None


def check_output(idx, test):
    msg0 = (
        "The output datasets should not be the same length.  There are "
        "15 reaches specified, but underlain by only 12 active cells. "
        "3 cell below the SFR/SFT reaches have idomain=0 and are "
        "therefore inactive."
    )
    msg1 = (
        "There seems to be misalignment in the flow/transport results "
        "when some of the GWF cells are inactive but overlain with SFR "
        "reaches"
    )
    msg2 = (
        "Calculated number of reaches with flow should be greater than "
        "the number of cells that receive an influx of mass from the "
        "losing reaches since 3 of the cells underlying the streams are "
        "inactive (idomain=0)."
    )

    sim = flopy.mf6.MFSimulation.load(sim_ws=test.workspace)
    gwf = sim.get_model()
    gwt = sim.get_model("gwt-" + cases[idx])
    sfr = gwf.get_package("SFR-1")
    sft = gwt.get_package("SFT")
    stage = sfr.output.stage().get_alldata().squeeze()
    sft_mass_loss = sft.output.budget().get_data(text="GWF")

    # The stream runs dry toward its bottom end, checking to make sure
    # that where the stream is dry, there is no simulated mass seepage
    # into the aquifer.  Where the stream is losing to the aquifer, checking
    # to ensure that there is some mass loss to the aquifer.
    for sp in np.arange(nper):
        assert len(sft_mass_loss[sp]) != len(stage[0]), msg0

        for i, mass_loss in enumerate(sft_mass_loss[sp]):
            gwf_cell_id = mass_loss[0]
            stg = stage[0][gwf_cell_id - 1]
            if stg > 0:
                assert mass_loss[-1] != 0.0, msg1
            else:
                assert np.isclose(mass_loss[-1], 0.0), msg1

        num_stg_vals = 0
        num_sftgwf_vals = 0
        for i, itm in enumerate(stage[sp]):
            if itm != 0:
                num_stg_vals += 1
        for i, itm in enumerate(sft_mass_loss[sp]):
            if not np.isclose(itm[-1], 0):
                num_sftgwf_vals += 1

        # Owing to 3 reaches being underlain by inactive cells and therefore
        # not able to receive mass from sft, the following values should not
        # be equal
        assert num_stg_vals > num_sftgwf_vals, msg2


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

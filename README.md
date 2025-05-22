
# MODFLOW 6: USGS Modular Hydrologic Model

This is the development repository for the USGS MODFLOW 6 Hydrologic Model. The  official USGS distribution is available at [USGS Release Page](https://water.usgs.gov/ogw/modflow/MODFLOW.html).

### Version 6.7.0.dev2 (preliminary)

[![GitHub release](https://img.shields.io/github/release/MODFLOW-ORG/modflow6.svg)](https://github.com/MODFLOW-ORG/modflow6/releases/latest)
[![MODFLOW 6 continuous integration](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/ci.yml/badge.svg)](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/ci.yml)
[![MODFLOW 6 documentation](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/docs.yml/badge.svg)](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/docs.yml)
[![MODFLOW 6 large models](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/large.yml/badge.svg)](https://github.com/MODFLOW-ORG/modflow6/actions/workflows/large.yml)
[![MODFLOW 6 nightly build](https://github.com/MODFLOW-ORG/modflow6-nightly-build/actions/workflows/nightly-build.yml/badge.svg)](https://github.com/MODFLOW-ORG/modflow6-nightly-build/actions/workflows/nightly-build.yml)

## Branches

This repository contains branches of ongoing MODFLOW 6 development.  The two main branches in this repository are:

* `master`: the state of the MODFLOW 6 repository corresponding to the last official USGS release
* `develop`: the current development version of the MODFLOW 6 program

The `develop` branch is under active and frequent updates by the MODFLOW development team and other interested contributors.  We follow a fork and pull request workflow and require that pull requests pass our test suite before they are considered a possible candidate to merge into `develop`. The `master` branch is only updated immediately prior to each new release.

This repository may contain other branches with various levels of development code; however, these branches may be merged into develop or deleted without notice.

## Nightly Builds

This repository's `develop` branch often contains bug fixes and new features that are not yet part of the official USGS release. Binaries for Linux, macOS and Windows are built from the `develop` branch and posted to the [`MODFLOW-ORG/modflow6-nightly-build` repository](https://github.com/MODFLOW-ORG/modflow6-nightly-build/releases) each night. The updated user guide `mf6io.pdf` is also included, as well as the `code.json` metadata file.

## Releases

Software distributions for the current and previous official USGS releases are available [here](https://github.com/MODFLOW-ORG/modflow6/releases).  The current release is also available [here](https://water.usgs.gov/water-resources/software/MODFLOW-6/) from the USGS.

## Examples
MODFLOW 6 has an extensive suite of example problems that are constructed using the python [FloPy](https://github.com/modflowpy/flopy) package.  These example problems are contained in a separate git repository located [here](https://github.com/MODFLOW-ORG/modflow6-examples).  These examples are included in the official [USGS MODFLOW 6 distribution](https://water.usgs.gov/water-resources/software/MODFLOW-6/), and they are also rendered into online [descriptions](https://modflow6-examples.readthedocs.io/en/master/examples.html) and [jupyter notebooks](https://modflow6-examples.readthedocs.io/en/master/notebook_examples.html).

## Contributing

First, follow [DEVELOPER.md](DEVELOPER.md) in order to set up your development environment to build and test MODFLOW 6.
The contribution workflow is described in [CONTRIBUTING.md](CONTRIBUTING.md).

## Continuous Integration

This repository contains an `./autotest` folder with python scripts for building and testing the MODFLOW 6 program and other workflow tasks.  The testing workflow relies heavily on several related repositories including:

* [modflowpy/pymake](https://github.com/modflowpy/pymake)
* [modflowpy/flopy](https://github.com/modflowpy/flopy)
* [MODFLOW-ORG/modflow6-examples](https://github.com/MODFLOW-ORG/modflow6-examples)
* [MODFLOW-ORG/modflow6-testmodels](https://github.com/MODFLOW-ORG/modflow6-testmodels)
* [MODFLOW-ORG/modflow6-largetestmodels](https://github.com/MODFLOW-ORG/modflow6-largetestmodels)
* [MODFLOW-ORG/executables](https://github.com/MODFLOW-ORG/executables)
* [MODFLOW-ORG/modflow-devtools](https://github.com/MODFLOW-ORG/modflow-devtools)
* [Deltares/xmipy](https://github.com/Deltares/xmipy)

## Code Documentation

Documentation for the code on the `develop` branch is automatically developed using Doxygen as part of the MODFLOW 6 Continuous Integration/Continuous Deployment process. Code documentation is available from the [GitHub Pages](https://MODFLOW-ORG.github.io/modflow6/) for this repository.

## Using Flopy with the `develop` Version of MODFLOW 6

All new MODFLOW 6 functionality in the `develop` branch can be used directly with flopy; however, the installed version of flopy may need to be updated in order to become aware of new MODFLOW 6 changes.  All MODFLOW 6 input is described using "definition files," which are located in the MODFLOW 6 repository under `./doc/mf6io/mf6ivar/dfn`.  These definition files may be updated on the `develop` branch as input changes to accommodate new options. The python script `./autotest/update_flopy.py` can be run to automatically update the installed version of flopy to work with the `develop` version of MODFLOW 6.  As a note of caution, be aware that running the `update_flopy.py` script will change the installed version of flopy, which could adversely affect existing MODFLOW 6 flopy scripts.

## Instructions for building definition files for new packages

Instructions for building definition files for new packages are summarized in [doc/mf6io/mf6ivar/readme.md](doc/mf6io/mf6ivar/readme.md).

## MODFLOW 6 Overview

MODFLOW is a popular open-source groundwater flow model distributed by the U.S. Geological Survey.  For 30 years, the MODFLOW program has been widely used by academic researchers, private consultants, and government scientists to accurately, reliably, and efficiently simulate groundwater flow.  With time, growing interest in surface and groundwater interactions, local refinement with nested and unstructured grids, karst groundwater flow, solute transport, and saltwater intrusion, has led to the development of numerous MODFLOW versions.  Although these MODFLOW versions are often based on the core version (presently MODFLOW-2005), there are often incompatibilities that restrict their use with one another.  In many cases, development of these alternative versions has been challenging due to the underlying MODFLOW structure, which was designed for the simulation with a single groundwater flow model using a rectilinear grid.

MODFLOW 6 is the latest core version of MODFLOW. It synthesizes many of the capabilities available in MODFLOW-2005, MODFLOW-NWT, and MODFLOW-LGR. MODFLOW 6 was built on a new object-oriented framework that allows new packages and models to be added, and allows any number of models to be run simultaneously in a single simulation.  Model may be coupled sequentially, such as for flow and transport, or the models may be tightly coupled at the matrix level, such as for multiple flow models. MODFLOW 6 presently contains two types of hydrologic models, the Groundwater Flow (GWF) Model and the Groundwater Transport (GWT) Model.

The Groundwater Flow (GWF) Model was the first model to be released in MODFLOW 6. It supports regular MODFLOW grids consisting of layers, rows, and columns, but it also supports more flexible grids that may conform to irregular boundaries or have increased resolution in areas of interest. The GWF Model consists of the original MODFLOW stress packages (CHD, WEL, DRN, RIV, GHB, RCH, and EVT) and four advanced stress packages (MAW, SFR, LAK, and UZF), which have been distilled from their predecessors to contain the most commonly used capabilities. MODFLOW 6 contains a new Water Mover (MVR) Package that can transfer water from provider packages to receiver packages. Providers can be many of the stress and advanced stress packages; receivers can be any of the advanced stress packages. This new capability makes it possible to route water between lakes and streams, route rejected infiltration into a nearby stream, or augment lakes using groundwater pumped from wells, for example. To modernize user interaction with the program, the MODFLOW 6 input structure was redesigned. Within package input files, information is divided into blocks, and informative keywords are used to label numeric data and activate options. This new input structure was designed to make it easier for users to adjust simulation options in an intuitive manner, reduce user input errors, and allow new capabilities to be added without causing problems with backward compatibility.

The GWT model for MODFLOW 6 simulates three-dimensional transport of a single solute species in flowing groundwater. The GWT Model solves the solute transport equation using numerical methods and a generalized CVFD approach, which can be used with regular MODFLOW grids or with unstructured grids. The GWT Model is designed to work with most of the new capabilities released with the GWF Model, including the Newton flow formulation, unstructured grids, advanced packages, and the movement of water between packages. The GWF and GWT Models operate simultaneously during a MODFLOW 6 simulation to represent coupled groundwater flow and solute transport. The GWT Model can also run separately from a GWF Model by reading the heads and flows saved by a previously run GWF Model. The GWT model is also capable of working with the flows from another groundwater flow model, as long as the flows from that model can be written in the correct form to flow and head files.


## How to Cite MODFLOW 6

### ***Report and Papers***

[Hughes, J.D., Langevin, C.D., and Banta, E.R., 2017, Documentation for the MODFLOW 6 framework: U.S. Geological Survey Techniques and Methods, book 6, chap. A57, 40 p., https://doi.org/10.3133/tm6A57.](https://doi.org/10.3133/tm6A57)

[Langevin, C.D., Hughes, J.D., Banta, E.R., Niswonger, R.G., Panday, Sorab, and Provost, A.M., 2017, Documentation for the MODFLOW 6 Groundwater Flow Model: U.S. Geological Survey Techniques and Methods, book 6, chap. A55, 197 p., https://doi.org/10.3133/tm6A55.](https://doi.org/10.3133/tm6A55)

[Provost, A.M., Langevin, C.D., and Hughes, J.D., 2017, Documentation for the "XT3D" option in the Node Property Flow (NPF) Package of MODFLOW 6: U.S. Geological Survey Techniques and Methods, book 6, chap. A56, 40 p., https://doi.org/10.3133/tm6A56.](https://doi.org/10.3133/tm6A56)

[Langevin, C.D., Panday, S, and Provost, A.M., 2020, Hydraulic-head formulation for density-dependent flow and transport: Groundwater, v. 58, no. 3, p. 349–362, https://doi.org/10.1111/gwat.12967.](https://doi.org/10.1111/gwat.12967)

[Morway, E.D., Langevin, C.D., and Hughes, J.D., 2021, Use of the MODFLOW 6 water mover package to represent natural and managed hydrologic connections: Groundwater, v. 59, no. 6, p. 913-924, https://doi.org/10.1111/gwat.13117.](https://doi.org/10.1111/gwat.13117)

[Hughes, J.D., Russcher, M.J., Langevin, C.D., Morway, E.D., and McDonald, R.R., 2022, The MODFLOW Application Programming Interface for simulation control and software interoperability: Environmental Modelling & Software, v. 148, 105257, https://doi.org/10.1016/j.envsoft.2021.105257.](https://doi.org/10.1016/j.envsoft.2021.105257)

[Langevin, C.D., Provost, A.M., Panday, Sorab, and Hughes, J.D., 2022, Documentation for the MODFLOW 6 Groundwater Transport (GWT) Model: U.S. Geological Survey Techniques and Methods, book 6, chap. A61, 56 p., https://doi.org/10.3133/tm6A61.](https://doi.org/10.3133/tm6A61)

[Hughes, J.D., Leake, S.A., Galloway, D.L., and White, J.T., 2022, Documentation for the Skeletal Storage, Compaction, and Subsidence (CSUB) Package of MODFLOW 6: U.S. Geological Survey Techniques and Methods, book 6, chap. A62, 57 p., https://doi.org/10.3133/tm6A62.](https://doi.org/10.3133/tm6A62)

[Langevin, C.D., Hughes, J.D., Provost, A.M., Russcher, M.J. and Panday, S., 2023, MODFLOW as a Configurable Multi-Model Hydrologic Simulator: Groundwater, v. 62, no. 1, p. 111-123, https://doi.org/10.1111/gwat.13351.](https://doi.org/10.1111/gwat.13351)

[Larsen, J.D., Langevin, C.D., Hughes, J.D. and Niswonger, R.G., 2024, An Agricultural Package for MODFLOW 6 Using the Application Programming Interface. Groundwater, v. 62, no. 1, p. 157-166, https://doi.org/10.1111/gwat.13367.](https://doi.org/10.1111/gwat.13367)

[Provost, A.M., Bardot, K., Langevin, C.D., McCallum, J.L., 2025, Accurate Simulation of Flow through Dipping Aquifers with MODFLOW 6 Using Enhanced Cell Connectivity. Groundwater, v. 63, no. 3, p. 399-408, https://doi.org/10.1111/gwat.13459.]( https://doi.org/10.1111/gwat.13459)

[Morway, E.D., Provost, A.M., Langevin, C.D., Hughes, J.D., Russcher, M.J., Chen, C.Y., Lin, Y.F.F., 2025, A New Groundwater Energy Transport Model for the MODFLOW Hydrologic Simulator. Groundwater, v. 63, no. 3, p. 409-421, https://doi.org/10.1111/gwat.13470.]( https://doi.org/10.1111/gwat.13470)

#### ***Software/Code***

The following is the general citation for the MODFLOW 6 software.

[Langevin, C.D., Hughes, J.D., Banta, E.R., Provost, A.M., Niswonger, R.G., and Panday, Sorab, 2017, MODFLOW 6 Modular Hydrologic Model: U.S. Geological Survey Software, https://doi.org/10.5066/F76Q1VQV](https://doi.org/10.5066/F76Q1VQV)

Citations for specific versions are included with the [releases](https://github.com/MODFLOW-ORG/modflow6/releases).


Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.


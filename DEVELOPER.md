# Developing MODFLOW 6

This document describes how to set up a development environment to modify, build and test MODFLOW 6. Details on how to contribute your code to the repository are found in the separate document [CONTRIBUTING.md](CONTRIBUTING.md). 

To build and test a parallel version of the program, first read the instructions below and then continue in [PARALLEL.md](PARALLEL.md).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Prerequisites](#prerequisites)
  - [Git](#git)
  - [Fortran compiler](#fortran-compiler)
    - [GNU Fortran](#gnu-fortran)
      - [Linux](#linux)
      - [macOS](#macos)
      - [Windows](#windows)
    - [Intel Fortran](#intel-fortran)
      - [Windows](#windows-1)
  - [Python](#python)
    - [Dependencies](#dependencies)
      - [`meson`](#meson)
      - [`fprettify`](#fprettify)
      - [`mfpymake`](#mfpymake)
      - [`flopy`](#flopy)
      - [`modflow-devtools`](#modflow-devtools)
  - [Optional tools](#optional-tools)
    - [GNU Make](#gnu-make)
    - [Visual Studio](#visual-studio)
    - [Doxygen & LaTeX](#doxygen--latex)
- [Installation](#installation)
- [Building](#building)
- [Testing](#testing)
  - [Configuring a test environment](#configuring-a-test-environment)
    - [Building development binaries](#building-development-binaries)
    - [Rebuilding and installing release binaries](#rebuilding-and-installing-release-binaries)
    - [Updating `flopy` plugins](#updating-flopy-plugins)
    - [External model repositories](#external-model-repositories)
    - [Installing external repos](#installing-external-repos)
      - [Test models](#test-models)
      - [Example models](#example-models)
  - [Running Tests](#running-tests)
    - [Selecting tests with markers](#selecting-tests-with-markers)
    - [External model tests](#external-model-tests)
    - [Writing tests](#writing-tests)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Prerequisites

Before you can build and test MODFLOW 6, you must install and configure the
following on your development machine:

- git
- Python3.8+
- a modern Fortran compiler

Some additional, optional tools are also discussed below.

### Git

[Git](https://git-scm.com) and/or the **GitHub app** (for [Mac](https://mac.github.com) or [Windows](https://windows.github.com)).
[GitHub's Guide to Installing Git](https://help.github.com/articles/set-up-git) is a good source of information.

### Fortran compiler

The GNU Fortran compiler `gfortran` or the Intel Fortran Classic compiler `ifort` can be used to compile MODFLOW 6.

**Note:** the next-generation Intel Fortran compiler `ifx` is not yet compatible with MODFLOW 6.

#### GNU Fortran

GNU Fortran can be installed on all three major platforms.

##### Linux

- fedora-based: `dnf install gcc-gfortran`
- debian-based: `apt install gfortran`

##### macOS

- [Homebrew](https://brew.sh/): `brew install gcc`
- [MacPorts](https://www.macports.org/): `sudo port install gcc10`

##### Windows

- Download the Minimalist GNU for Windows (MinGW) installer from Source Forge:
  https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe
- Run the installer. Make sure to change `Architecture` to `x86_64`. Leave the
  other settings on default.
- Find the `mingw64/bin` directory in the installation and add it
  to your PATH. Find `Edit the system environment variables` in your Windows
  Start Screen. Click the `Environmental Variables` button and double-click the
  `Path` variable in the User Variables (the top table). Click the `New` button
  and enter the location of the `mingw64/bin` directory.

#### Intel Fortran

Intel Fortran can also be used to compile MODFLOW 6 and associated utilities. The `ifort` compiler is available in the [Intel oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html). An installer is bundled with the download.

A number of environment variables must be set before using Intel Fortran. General information can be found [here](https://www.intel.com/content/www/us/en/develop/documentation/oneapi-programming-guide/top/oneapi-development-environment-setup.html), with specific instructions to configure a shell session for `ifort` [here](https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/compiler-setup/use-the-command-line/specifying-the-location-of-compiler-components.html).

##### Windows

On Windows, [Visual Studio](https://visualstudio.microsoft.com) and a number of libraries must be installed for `ifort` to work. The required libraries can be installed by ticking the "Desktop Development with C++" checkbox in the Visual Studio Installer's Workloads tab. 

**Note:** Invoking the `setvars.bat` scripts from a Powershell session will *not* put `ifort` on the path, since [batch script environments are local to their process](https://stackoverflow.com/a/49028002/6514033). To relaunch PowerShell with oneAPI variables configured:

```
cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" && "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" && powershell'
```

### Python

Python 3.8+ is required to run MODFLOW 6 tests. A Conda distribution (e.g. [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/individual) is recommended. Python dependencies are specified in `environment.yml`. To create an environment, run from the project root:

```
conda env create -f environment.yml
```

To update an existing environment:

```shell
conda env update -f environment.yml
```

#### Dependencies

This project depends critically on a few Python packages for building, linting and testing tasks:

- `meson`
- `fprettify`
- `pymake`
- `flopy`

These are each described briefly below. The Conda `environment.yml` contains a number of other dependencies also required for various development tasks, but they are not described in detail here.

##### `meson`

[Meson](https://mesonbuild.com/index.html) is recommended for building MODFLOW 6 and is included in `environment.yml`. It can also be [installed independently](https://mesonbuild.com/Getting-meson.html) &mdash; note that if you do so you will need to manually add the executable to the [PATH](https://en.wikipedia.org/wiki/PATH_(variable)).

##### `fprettify`

[`fprettify`](https://github.com/pseewald/fprettify) can be used to format Fortran source code and in combination with the [MODFLOW 6 fprettify configuration](https://github.com/MODFLOW-USGS/modflow6/blob/develop/distribution/.fprettify.yaml) establishes a contribution standard for properly formatted MODFLOW 6 Fortran source. This tool can be installed with `pip` or `conda` and used from the command line or integrated with a [VSCode](https://github.com/MODFLOW-USGS/modflow6/blob/develop/.vscode/README.md) or Visual Studio development environment. The `fprettify` package is included in the Conda environment in `environment.yml`. See [contribution guidelines](https://github.com/MODFLOW-USGS/modflow6/blob/develop/CONTRIBUTING.md) for additional information.

##### `mfpymake`

The `mfpymake` package can build MODFLOW 6 and related programs and artifacts (e.g. makefiles), and is used in particular by the `distribution/build_makefiles.py` script. `mfpymake` is included in the Conda environment in `environment.yml`. To install separately, follow the instructions as explained on the README of the [repository](https://github.com/modflowpy/pymake). The README also demonstrates basic usage.

##### `flopy`

[`flopy`](https://github.com/modflowpy/flopy) is used throughout MODFLOW 6 tests to create, run and post-process models.

Like MODFLOW 6, `flopy` is modular &mdash; for each MODFLOW 6 package there is generally a corresponding `flopy` plugin. Plugins are generated dynamically from DFN files stored in this repository under `doc/mf6io/mf6ivar/dfn`.

##### `modflow-devtools`

The tests use a set of shared fixtures and utilities provided by the [`modflow-devtools`](https://github/com/MODFLOW-USGS/modflow-devtools) package. This package is included in the Conda environment in `environment.yml`.

### Optional tools

Some other tools are useful but not required to develop MODFLOW 6.

#### GNU Make

This repository provides makefiles, generated by `mfpymake`, which can be used to build MODFLOW 6 with [GNU Make](https://www.gnu.org/software/make/). For further instructions we refer to the [GNU Make Manual](https://www.gnu.org/software/make/manual/).

#### Visual Studio

Visual Studio installers can be downloaded from the [official website](https://visualstudio.microsoft.com/). MODFLOW 6 solution files can be found in the `msvs` folder.

#### Doxygen & LaTeX

[Doxygen](https://www.doxygen.nl/index.html) is used to generate the [MODFLOW 6 source code documentation](https://modflow-usgs.github.io/modflow6/). [Graphviz](https://graphviz.org/) is used by doxygen to produce source code diagrams. [LaTeX](https://www.latex-project.org/) is used to generate the MODFLOW 6 release notes and Input/Output documents (docs/mf6io/mf6io.nightlybuild).

These programs can be installed from various sources, including by conda, macports, or from individual sources such as https://www.tug.org/. Details about USGS LaTeX libraries can be seen in addition to linux installs in the CI workflow for the docs (`.github/workflows/ci-docs.yml`).

## Installation

Fork and clone the MODFLOW 6 repository:

1. Login to your GitHub account or create one by following the instructions given [here](https://github.com/signup/free).
2. [Fork](http://help.github.com/forking) the [main MODFLOW 6](https://github.com/MODFLOW-USGS/modflow6).
3. Clone your fork of the MODFLOW 6 repository and create an `upstream` remote pointing back to your fork.

```shell
# Clone your GitHub repository:
git clone git@github.com:<github username>/modflow6.git

# Go to the MODFLOW 6 directory:
cd modflow6

# Add the main MODFLOW 6 repository as an upstream remote to your repository:
git remote add upstream https://github.com/MODFLOW-USGS/modflow6.git
```

## Building

Meson is the recommended build tool for MODFLOW 6. [Meson](https://mesonbuild.com/Getting-meson.html) must be installed and on your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)). Creating and activating the Conda environment `environment.yml` should be sufficient for this.

Meson build configuration files are provided for MODFLOW 6 as well as `zbud6` and `mf5to6` utility programs:

- `meson.build`
- `utils/zonebudget/meson.build`
- `utils/mf5to6/meson.build`

To build MODFLOW 6, first configure the build directory. By default Meson uses compiler flags for a release build. To create a debug build, add `-Doptimization=0` to the following `setup` command.

```shell
# bash (linux and macOS)
meson setup builddir --prefix=$(pwd) --libdir=bin

# cmd (windows)
meson setup builddir --prefix=%CD% --libdir=bin
```

Compile MODFLOW 6 by executing:

```shell
meson compile -C builddir
```

In order to run the tests the binaries have to be installed:

```shell
meson install -C builddir
```

The binaries can then be found in the `bin` folder. `meson install` also triggers a compilation if necessary, so executing `meson install` is enough to get up-to-date binaries in the `bin` folder.

**Note:** If using Visual Studio Code, you can use tasks as described [here](.vscode/README.md) to automate the above.

## Testing

MODFLOW 6 tests are driven with [`pytest`](https://docs.pytest.org/en/7.1.x/), with the help of plugins like `pytest-xdist` and `pytest-cases`. Testing dependencies are included in the Conda environment `environment.yml`.

**Note:** the entire test suite should pass before a pull request is submitted. Tests run in GitHub Actions CI and a PR can only be merged with passing tests. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for more information.

### Configuring a test environment

A few tasks must be completed before running tests:

- build local MODFLOW 6 development version
- rebuild the last MODFLOW 6 release
- install additional executables
- update FloPy packages and plugins
- clone MODFLOW 6 test model and example repositories

Tests expect binaries to live in the `bin` directory relative to the project root, as configured above in the `meson` commands. Binaries are organized as follows:

- local development binaries in the top-level `bin` folder
- binaries rebuilt in development mode from the latest release in `bin/rebuilt`
- related programs installed from the [executables distribution](https://github.com/MODFLOW-USGS/executables/releases) live in `bin/downloaded`

Tests must be run from the `autotest` folder.

#### Building development binaries

Before running tests, the local development version of MODFLOW 6 must be built with `meson` as described above. The `autotest/build_exes.py` script is provided as a shortcut to easily rebuild local binaries. The script can be run from the project root with:

```shell
python autotest/build_exes.py
```

Alternatively, it can be run from the `autotest` directory with `pytest`:

```shell
pytest build_exes.py
```

By default, binaries will be placed in the `bin` directory relative to the project root, as in the `meson` commands described above. To change the location of the binaries, use the `--path` option.

#### Rebuilding and installing release binaries

Tests require the latest official MODFLOW 6 release to be compiled in develop mode with the same Fortran compiler as the development version. A number of binaries distributed from the [executables repo](https://github.com/MODFLOW-USGS/executables) must also be installed. The script `autotest/get_exes.py` does both of these things. It can be run from the project root with:

```shell
python autotest/get_exes.py
```

Alternatively, with `pytest` from the `autotest` directory:

```shell
pytest get_exes.py
```

By default, binaries will be placed in the `bin` directory relative to the project root, as in the `meson` commands described above. Nested `bin/downloaded` and `bin/rebuilt` directories are created to contain the rebuilt last release and the downloaded executables, respectively. To change the location of the binaries, use the `--path` option.

#### Updating `flopy` plugins

Plugins should be regenerated from DFN files before running tests for the first time or after definition files change. This can be done with the `autotest/update_flopy.py` script, which wipes and regenerates plugin classes for the `flopy` installed in the Python environment.

**Note:** if you've installed a local version of `flopy` from source, running this script can overwrite files in your repository.

There is a single optional argument, the path to the folder containing definition files. By default DFN files are assumed to live in `doc/mf6io/mf6ivar/dfn`, making the following identical:

```shell
python autotest/update_flopy.py
python autotest/update_flopy.py doc/mf6io/mf6ivar/dfn
```

#### External model repositories

Some autotests load example models from external repositories:

- [`MODFLOW-USGS/modflow6-testmodels`](https://github.com/MODFLOW-USGS/modflow6-testmodels)
- [`MODFLOW-USGS/modflow6-largetestmodels`](https://github.com/MODFLOW-USGS/modflow6-largetestmodels)
- [`MODFLOW-USGS/modflow6-examples`](https://github.com/MODFLOW-USGS/modflow6-examples)

#### Installing external repos

By default, the tests expect these repositories side-by-side with (i.e. in the same parent directory as) the `modflow6` repository. If the repos are somewhere else, you can set the `REPOS_PATH` environment variable to point to their parent directory. If external model repositories are not found, tests requiring them will be skipped.

**Note:** a convenient way to persist environment variables needed for tests is to store them in a `.env` file in the `autotest` folder. Each variable should be defined on a separate line, with format `KEY=VALUE`. The `pytest-dotenv` plugin will then automatically load any variables found in this file into the test process' environment.

##### Test models

The test model repos can simply be cloned &mdash; ideally, into the parent directory of the `modflow6` repository, so that repositories live side-by-side:

```shell
git clone MODFLOW-USGS/modflow6-testmodels
git clone MODFLOW-USGS/modflow6-largetestmodels
```

##### Example models

First clone the example models repo:

```shell
git clone MODFLOW-USGS/modflow6-examples
```

The example models require some setup after cloning. Some extra Python dependencies are required to build the examples: 

```shell
cd modflow6-examples/etc
pip install -r requirements.pip.txt
```

Then, still from the `etc` folder, run:

```shell
python ci_build_files.py
```

This will build the examples for subsequent use by the tests.

### Running Tests

Tests are driven by `pytest` and must be run from the `autotest` folder. To run tests in a particular file, showing verbose output, use:

```shell
pytest -v <file>
```

Tests can be run in parallel with the `-n` option, which accepts an integer argument for the number of parallel processes. If the value `auto` is provided, `pytest-xdist` will use one worker per available processor.

```shell
pytest -v -n auto
```

#### Selecting tests with markers

Markers can be used to select subsets of tests. Markers provided in `pytest.ini` include:

- `slow`: tests that take longer than a few seconds to complete
- `repo`: tests that require external model repositories
- `large`: tests using large models (from the `modflow6-examples` and `modflow6-largetestmodels` repos)
- `regression`: tests comparing results from multiple versions

Markers can be used with the `-m <marker>` option, and can be applied in boolean combinations with `and`, `or` and `not`. For instance, to run fast tests in parallel, excluding regression tests:

```shell
pytest -v -n auto -m "not slow and not regression"
```

The `--smoke` (short `-S`) flag, provided by `modflow-devtools` is an alias for the above:

```shell
pytest -v -n auto -S
```

[Smoke testing](https://modflow-devtools.readthedocs.io/en/latest/md/markers.html#smoke-testing) is a form of integration testing which aims to test a decent fraction of the codebase quickly enough to run often during development.

#### External model tests

Tests using models from external repositories can be selected with the `repo` marker:

```shell
pytest -v -n auto -m "repo"
```

The `large` marker is a subset of the `repo` marker. To test models excluded from commit-triggered CI and only run on GitHub Actions nightly:

```shell
pytest -v -n auto -m "large"
```

Test scripts for external model repositories can also be run independently:

```shell
# MODFLOW 6 test models
pytest -v -n auto test_z01_testmodels_mf6.py

# MODFLOW 5 to 6 conversion test models
pytest -v -n auto test_z02_testmodels_mf5to6.py

# models from modflow6-examples repo
pytest -v -n auto test_z03_examples.py

# models from modflow6-largetestmodels repo
pytest -v -n auto test_z03_largetestmodels.py
```

Tests load external models from fixtures provided by `modflow-devtools`. External model tests can be selected by model or simulation name, or by packages used. See the [`modflow-devtools` documentation](https://modflow-devtools.readthedocs.io/en/latest/md/fixtures.html#filtering) for usage examples. Note that filtering options only apply to tests using external models, and will not filter tests defining models in code &mdash; for that, the `pytest` built-in `-k` option may be used.

#### Writing tests

Tests should ideally follow a few conventions for easier maintenance:

- Use temporary directory fixtures. Tests which write to disk should use `pytest`'s built-in `tmp_path` fixtures or one of the [keepable temporary directory fixtures from `modflow-devtools`](https://modflow-devtools.readthedocs.io/en/latest/md/fixtures.html#keepable-temporary-directories). This prevents tests from polluting one another's state.

- Use markers for convenient (de-)selection:
  - `@pytest.mark.slow` if the test doesn't complete in a few seconds (this preserves the ability to quickly [`--smoke` test](https://modflow-devtools.readthedocs.io/en/latest/md/markers.html#smoke-testing)
  - `@pytest.mark.repo` if the test relies on external model repositories
  - `@pytest.mark.regression` if the test compares results from different versions

**Note:** If all three external model repositories are not installed as described above, some tests will be skipped. The full test suite includes >750 cases. All must pass before changes can be merged into this repository.

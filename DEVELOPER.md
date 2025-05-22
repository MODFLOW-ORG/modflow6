# Developing MODFLOW 6

This document describes how to set up a development environment to modify, build and test MODFLOW 6. Details on how to contribute your code to the repository are found in the separate document [CONTRIBUTING.md](./CONTRIBUTING.md). 

To build and test an extended version of the program, first read the instructions below and then continue in [EXTENDED.md](./EXTENDED.md).

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Prerequisites](#prerequisites)
  - [Git](#git)
  - [Python](#python)
  - [Fortran compiler](#fortran-compiler)
    - [GNU Fortran](#gnu-fortran)
    - [Intel Fortran](#intel-fortran)
      - [Windows](#windows)
    - [Compiler compatibility](#compiler-compatibility)
      - [Compile](#compile)
      - [Test](#test)
  - [Optional tools](#optional-tools)
- [Get the MODFLOW 6 repository](#get-the-modflow-6-repository)
- [Install the python environment](#install-the-python-environment)
  - [Python dependencies](#python-dependencies)
    - [`meson`](#meson)
    - [`codespell`](#codespell)
    - [`fprettify`](#fprettify)
    - [`ruff`](#ruff)
    - [`mfpymake`](#mfpymake)
    - [`flopy`](#flopy)
    - [`modflow-devtools`](#modflow-devtools)
- [Building](#building)
- [Formatting](#formatting)
  - [Spell checking](#spell-checking)
  - [Fortran formatting](#fortran-formatting)
  - [Python formatting](#python-formatting)
  - [Python linting](#python-linting)
- [Testing](#testing)
  - [Configuring a test environment](#configuring-a-test-environment)
    - [Configuring unit tests](#configuring-unit-tests)
    - [Configuring integration tests](#configuring-integration-tests)
      - [Rebuilding release binaries](#rebuilding-release-binaries)
      - [Updating FloPy packages](#updating-flopy-packages)
      - [Updating Fortran definitions](#updating-fortran-definitions)
      - [Installing external models](#installing-external-models)
  - [Running tests](#running-tests)
    - [Running unit tests](#running-unit-tests)
    - [Running integration tests](#running-integration-tests)
      - [Selecting tests with markers](#selecting-tests-with-markers)
  - [Writing tests](#writing-tests)
    - [Writing unit tests](#writing-unit-tests)
    - [Writing integration tests](#writing-integration-tests)
      - [Test framework](#test-framework)
- [Generating makefiles](#generating-makefiles)
  - [Updating extra and excluded files](#updating-extra-and-excluded-files)
  - [Testing makefiles](#testing-makefiles)
  - [Installing `make` on Windows](#installing-make-on-windows)
    - [Using Conda from Git Bash](#using-conda-from-git-bash)
- [Branching model](#branching-model)
  - [Overview](#overview)
  - [Managing long-lived branches](#managing-long-lived-branches)
    - [Backup](#backup)
    - [Squash](#squash)
    - [Rebase](#rebase)
    - [Cleanup](#cleanup)
- [Deprecation policy](#deprecation-policy)
  - [Finding deprecations](#finding-deprecations)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Prerequisites

Before you can build and test MODFLOW 6, you must install and configure the following on your development machine:

- git
- Python3.10+
- a modern Fortran compiler

Some additional, optional tools are also discussed below.

### Git

[Git](https://git-scm.com) and/or the **GitHub app** (for [Mac](https://mac.github.com) or [Windows](https://windows.github.com)).
[GitHub's Guide to Installing Git](https://help.github.com/articles/set-up-git) is a good source of information.

Optionally, the [`git blame`](https://git-scm.com/docs/git-blame) tool can be configured to work locally using:

```shell
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

### Python

Python 3.10+ is required to run MODFLOW 6 tests and in some cases to build MODFLOW 6. Information on installing the python environment is given in the [Installing Python environment](#install-the-python-environment) section. The MODFLOW 6 python environment should be installed after [locally cloning the repository](#get-the-modflow-6-repository).

### Fortran compiler

GNU Fortran or Intel Fortran compilers can be used to build MODFLOW 6. It may be possible to build MODFLOW 6 with other compilers, but this cannot be guaranteed.

#### GNU Fortran

GNU Fortran can be installed on all three major platforms.

*Linux*

- Fedora-based: `dnf install gcc-gfortran`
- Debian-based: `apt install gfortran`

*macOS*

- [Homebrew](https://brew.sh/): `brew install gcc@13`
- [MacPorts](https://www.macports.org/): `sudo port install gcc13`

**Note:** Xcode 15 includes a new linker implementation which breaks GNU Fortran compatibility. A workaround is to set `LDFLAGS` to use the classic linker, for instance:

```shell
export LDFLAGS="$LDFLAGS -Wl,-ld_classic"
```

See [this ticket](https://github.com/mesonbuild/meson/issues/12282) on the Meson repository for more information.

*Windows*

[Minimalist GNU for Windows](https://www.mingw-w64.org/) is the recommended way to obtain the GCC toolchain on Windows. Several MinGW distributions are available.

To install with Chocolatey: `choco install mingw`

To install from SourceForge:

- Download the MinGW installer:
  https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win32/Personal%20Builds/mingw-builds/installer/mingw-w64-install.exe
- Run the installer. Make sure to change `Architecture` to `x86_64`. Leave the
  other settings on default.
- Find the `mingw64/bin` directory in the installation and add it
  to your PATH. Find `Edit the system environment variables` in your Windows
  Start Screen. Click the `Environmental Variables` button and double-click the
  `Path` variable in the User Variables (the top table). Click the `New` button
  and enter the location of the `mingw64/bin` directory.

Binaries may also be downloaded and installed from the [releases here](https://github.com/brechtsanders/winlibs_mingw/releases).

**Note:** the MinGW distribution [available on conda-forge](https://anaconda.org/conda-forge/m2w64-toolchain_win-64) includes an outdated version of GCC and is not compatible with MODFLOW 6.

#### Intel Fortran

Intel Fortran can also be used to compile MODFLOW 6 and associated utilities. The `ifort` and `ifx` compilers are available in the [Intel oneAPI HPC Toolkit](https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html).

A number of environment variables must be set before using Intel Fortran. General information can be found [here](https://www.intel.com/content/www/us/en/develop/documentation/oneapi-programming-guide/top/oneapi-development-environment-setup.html), with specific instructions to configure a shell session for `ifort` [here](https://www.intel.com/content/www/us/en/develop/documentation/fortran-compiler-oneapi-dev-guide-and-reference/top/compiler-setup/use-the-command-line/specifying-the-location-of-compiler-components.html).

While the current development version of MODFLOW 6 is broadly compatible with `ifort`, `ifx` compatibility is still limited on Ubuntu and Windows, and `ifx` is not supported on macOS.

##### Windows

On Windows, [Visual Studio](https://visualstudio.microsoft.com) and a number of libraries must be installed for `ifort` and `ifx` to work. The required libraries can be installed by ticking the "Desktop Development with C++" checkbox in the Visual Studio Installer's Workloads tab. 

**Note:** Invoking the `setvars.bat` scripts from a Powershell session will *not* put `ifort` or `ifx` on the path, since [batch script environments are local to their process](https://stackoverflow.com/a/49028002/6514033). To relaunch PowerShell with oneAPI variables configured:

```
cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" && "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\env\vars.bat" && powershell'
```

#### Compiler compatibility

The following tables are automatically generated by [a CI workflow](.github/workflows/compilers.yml).

##### Compile

<!-- compile compat starts -->
| runner       | gcc 10   | gcc 11   | gcc 12   | gcc 13   | gcc 7   | gcc 8   | gcc 9   | intel-classic 2021.1   | intel-classic 2021.10   | intel-classic 2021.2   | intel-classic 2021.3   | intel-classic 2021.4   | intel-classic 2021.5   | intel-classic 2021.6   | intel-classic 2021.7   | intel-classic 2021.8   | intel-classic 2021.9   |   intel 2021.1 |   intel 2021.2 |   intel 2021.4 |   intel 2022.0 |   intel 2022.1 | intel 2022.2.1   | intel 2022.2   |   intel 2023.0 |   intel 2023.1 | intel 2023.2   |
|:-------------|:----------------|:----------------|:----------------|:----------------|:---------------|:---------------|:---------------|:------------------------------|:-------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|----------------------:|----------------------:|----------------------:|----------------------:|----------------------:|:------------------------|:----------------------|----------------------:|----------------------:|:----------------------|
| macos-11     | &check;         | &check;         | &check;         | &check;         | &check;        | &check;        | &check;        | &check;                       | &check;                        | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| macos-12     | &check;         | &check;         | &check;         | &check;         | &check;        | &check;        | &check;        | &check;                       | &check;                        | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| ubuntu-20.04 | &check;         | &check;         |              |              | &check;        | &check;        | &check;        | &check;                       | &check;                        | &check;                       |                            | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    | &check;                 | &check;               |                    |                    | &check;               |
| ubuntu-22.04 | &check;         | &check;         | &check;         | &check;         |             |             | &check;        | &check;                       | &check;                        | &check;                       |                            | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    | &check;                 | &check;               |                    |                    | &check;               |
| windows-2019 | &check;         | &check;         | &check;         | &check;         |             |             | &check;        |                            | &check;                        |                            |                            |                            |                            | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    |                      | &check;               |                    |                    | &check;               |
| windows-2022 | &check;         | &check;         | &check;         | &check;         |             |             | &check;        |                            | &check;                        |                            |                            |                            |                            | &check;                       | &check;                       | &check;                       | &check;                       |                    |                    |                    |                    |                    |                      | &check;               |                    |                    | &check;               |
<!-- compile compat ends -->

##### Test

<!-- test compat starts -->
| runner       | gcc 10   | gcc 11   | gcc 12   | gcc 13   | gcc 7   | gcc 8   | gcc 9   | intel-classic 2021.1   |   intel-classic 2021.10 | intel-classic 2021.2   | intel-classic 2021.3   | intel-classic 2021.4   | intel-classic 2021.5   | intel-classic 2021.6   | intel-classic 2021.7   |   intel-classic 2021.8 |   intel-classic 2021.9 |   intel 2021.1 |   intel 2021.2 |   intel 2021.4 |   intel 2022.0 |   intel 2022.1 |   intel 2022.2.1 |   intel 2022.2 |   intel 2023.0 |   intel 2023.1 |   intel 2023.2 |
|:-------------|:----------------|:----------------|:----------------|:----------------|:---------------|:---------------|:---------------|:------------------------------|-------------------------------:|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|:------------------------------|------------------------------:|------------------------------:|----------------------:|----------------------:|----------------------:|----------------------:|----------------------:|------------------------:|----------------------:|----------------------:|----------------------:|----------------------:|
| macos-11     | &check;         | &check;         | &check;         | &check;         | &check;        | &check;        | &check;        | &check;                       |                             | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| macos-12     | &check;         | &check;         | &check;         | &check;         |             |             |             | &check;                       |                             | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| ubuntu-20.04 | &check;         | &check;         |              |              | &check;        | &check;        | &check;        | &check;                       |                             | &check;                       |                            | &check;                       | &check;                       | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| ubuntu-22.04 | &check;         | &check;         | &check;         | &check;         |             |             | &check;        | &check;                       |                             | &check;                       |                            | &check;                       | &check;                       | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| windows-2019 |              |              |              | &check;         |             |             |             |                            |                             |                            |                            |                            |                            | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
| windows-2022 | &check;         | &check;         | &check;         | &check;         |             |             | &check;        |                            |                             |                            |                            |                            |                            | &check;                       | &check;                       |                            |                            |                    |                    |                    |                    |                    |                      |                    |                    |                    |                    |
<!-- test compat ends -->

### Optional tools

Some other tools are useful but not required to develop MODFLOW 6.

*GNU Make*

This repository provides makefiles, generated by `mfpymake`, which can be used to build MODFLOW 6 with [GNU Make](https://www.gnu.org/software/make/). For further instructions we refer to the [GNU Make Manual](https://www.gnu.org/software/make/manual/).

*Visual Studio*

Visual Studio installers can be downloaded from the [official website](https://visualstudio.microsoft.com/). MODFLOW 6 solution files can be found in the `msvs` folder.

*Doxygen & LaTeX*

[Doxygen](https://www.doxygen.nl/index.html) is used to generate the [MODFLOW 6 source code documentation](https://MODFLOW-ORG.github.io/modflow6/). [Graphviz](https://graphviz.org/) is used by doxygen to produce source code diagrams. [LaTeX](https://www.latex-project.org/) is used to generate the MODFLOW 6 release notes and Input/Output documents.

These programs can be installed from various sources, including by conda, macports, or from individual sources such as https://www.tug.org/. Details about USGS LaTeX libraries can be seen in addition to linux installs in the CI workflow for the docs (`.github/workflows/ci-docs.yml`).


## Get the MODFLOW 6 repository

Fork and clone the MODFLOW 6 repository:

1. Login to your GitHub account or create one by following the instructions given [here](https://github.com/signup/free).
2. [Fork](http://help.github.com/forking) the [main MODFLOW 6](https://github.com/MODFLOW-ORG/modflow6).
3. Clone your fork of the MODFLOW 6 repository and create an `upstream` remote pointing back to your fork.

After forking the MODFLOW 6 repository on GitHub.

1. Clone your fork of the GitHub repository to your computer.

```shell
  git clone git@github.com:<github username>/modflow6.git
```

2. Go to the MODFLOW 6 directory.

```shell
cd modflow6
```

3. Add the main MODFLOW 6 repository as an upstream remote to your repository.

```shell
git remote add upstream https://github.com/MODFLOW-ORG/modflow6.git
```

## Install the python environment

Python 3.10+ is required to run MODFLOW 6 tests and in some cases to build MODFLOW 6. Miniforge is the recommended python distribution if you do not have an existing Conda or Mamba based python distribution.

The [environment file for MODFLOW 6](./environment.yml) includes all of the required [python dependencies](#python-dependencies). Install the `modflow6` environment using the Conda `environment.yml` file in the repository. 

1. Open a terminal (command prompt) in the root directory of the repository.
2. Use either Mamba or Conda to install the `modflow6` environment.

```shell
mamba env create -f environment.yml 
```

```shell
conda env create -f environment.yml
```

Python can also be installed via Pixi. Pixi is currently being used to install python on GitHub Actions continuous integration/continuous development (CI/CD) virtual machines. In the future, Pixi may be the preferred approach for installing python for MODFLOW 6. As a result it is recommended for developers to also install the Pixi python environment, which can coexist with the Mamba/Conda python installation and `modflow6` environment. 

Pixi installation docs can be found [here](https://pixi.sh). After installing `pixi`, to set up an environment with all development dependencies, in the root directory of the MODFLOW 6 repository run:

```shell
pixi run install
```

### Python dependencies

This project depends critically on a few Python packages for building, linting, spell checking, and testing tasks:

- `meson`
- `codespell`
- `fprettify`
- `ruff`
- `mfpymake`
- `flopy`
- `modflow-devtools`

These are each described briefly below. These and a number of other dependencies are build-, test-, or release-time dependencies are included the Pixi environment `pixi.toml` as well as the Conda `environment.yml` file in this repository.

#### `meson`

[Meson](https://mesonbuild.com/index.html) is the recommended build system for MODFLOW 6.

#### `codespell`

[`codespell`](https://github.com/codespell-project/codespell) was designed primarily for checking misspelled words in source code, but can be used with other text files as well. This tool can be used from the command line or integrated with a [VSCode](.vscode/README.md). See [Spell check guidelines](#spell-checking) for additional information.

#### `fprettify`

[`fprettify`](https://github.com/pseewald/fprettify) can be used to format Fortran source code and in combination with the [MODFLOW 6 fprettify configuration](.fprettify.yaml) establishes a contribution standard for properly formatted MODFLOW 6 Fortran source. This tool can be used from the command line or integrated with a [VSCode](.vscode/README.md) or Visual Studio development environment. See [Fortran formatting guidelines](#fortran-formatting) for additional information.

#### `ruff`

[`ruff`](https://docs.astral.sh/ruff/) can be used to format and lint python code and scripts (for example, autotest scripts) and in combination with the [MODFLOW 6 ruff configuration](.github/common/ruff.toml) establishes a contribution standard for properly formatted python code and scripts. This tool can be used from the command line or integrated with a [VSCode](.vscode/README.md).  See [python formatting guidelines](#python-formatting) and [python linting guidelines](#python-linting) for additional information.

#### `mfpymake`

The `mfpymake` package can build MODFLOW 6 and related programs and artifacts (e.g. makefiles), and is used in particular by the `distribution/build_makefiles.py` script.

#### `flopy`

[`flopy`](https://github.com/modflowpy/flopy) is used throughout MODFLOW 6 tests to create, run and post-process models.

Like MODFLOW 6, `flopy` is modular &mdash; for each MODFLOW 6 package there is generally a corresponding `flopy` package. Packages are generated dynamically from DFN files stored in this repository under `doc/mf6io/mf6ivar/dfn`.

#### `modflow-devtools`

The tests use a set of shared fixtures and utilities provided by the [`modflow-devtools`](https://github.com/MODFLOW-ORG/modflow-devtools) package.

## Building

Meson is the recommended build tool for MODFLOW 6. [Meson](https://mesonbuild.com/Getting-meson.html) must be installed and on your [PATH](https://en.wikipedia.org/wiki/PATH_(variable)). Creating and activating the provided Pixi or Conda environment should be sufficient for this.

### MODFLOW 6 and ZONEBUDGET
Meson build configuration files are provided for MODFLOW 6 and the ZONEBUDGET utility program, and for Fortran unit tests (see [Testing](#testing) section below).

- `meson.build`
- `utils/zonebudget/meson.build`
- `autotest/meson.build`

Building MODFLOW 6 requires two steps:

- configure the build directory
- build the project

To configure the build directory for a debug version:

```shell
meson setup --prefix=$(pwd) --libdir=bin builddir -Ddebug=true
```
Or to configure the build directory for an optimized release version:

```shell
meson setup --prefix=$(pwd) --libdir=bin builddir
```

or using pixi to setup the build directory:

```shell
pixi run setup builddir
```

Debug versions can be built using pixi by adding `-Ddebug=true` at the end of the pixi command. Other meson commands (for example, `-Dextended=true`, `--wipe`, _etc._) added to the pixi command are passed through to Meson.

Substitute `%CD%` as necessary on Windows.

To build MODFLOW 6 and install binaries to `<project root>/bin/`:

```shell
meson install -C builddir
```

or using pixi:

```shell
pixi run build builddir
```

### MODFLOW 2005 to 6 converter
Meson build configuration files are provided for the MODFLOW 2005 to 6 converter utility program.

- `utils/mf5to6/meson.build`
- `utils/mf5to6/src/meson.build`

Building MODFLOW 2005 to 6 converter program requires two steps:

- configure the build directory
- build the project

To configure the build directory for a debug version from the `<project root>/utils/mf5to6` directory:

```shell
meson setup --prefix=$(pwd)/../../  builddir -Ddebug=true
```
Or to configure the build directory for an optimized release version from the `<project root>/utils/mf5to6` directory:

```shell
meson setup --prefix=$(pwd)/../../ builddir
```

or using pixi to setup the build directory from the `<project root>` directory:

```shell
pixi run setup-mf5to6 builddir
```

Debug versions can be built using pixi by adding `-Ddebug=true` at the end of the pixi command. Other meson commands (for example, `--wipe`, _etc._) added to the pixi command are passed through to Meson.

Substitute `%CD%` as necessary on Windows.

To build MODFLOW 6 and install binaries to `<project root>/bin/` from the `<project root>/utils/mf5to6` directory:

```shell
meson install -C builddir
```

or using pixi from the `<project root>` directory:

```shell
pixi run build-mf5to6 builddir
```

**Note:** If using Visual Studio Code, you can use tasks as described [here](.vscode/README.md) to automate the above.

## Formatting

### Spell checking

Fortran source files, python files, definition files, markdown, and LaTeX files can be checked with [codespell](https://github.com/codespell-project/codespell). codespell was designed primarily for checking misspelled words in source code, but it can be used with other text files as well. The `codespell` package is included in the Conda `environment.yml` and the Pixi `pixi.toml` files and can be run directly, via Pixi, or via [VSCode](.vscode/README.md) tasks.

To check whether the repository's Fortran source files, python files, definition files, markdown, and LaTeX files have any spelling errors without making any changes:

```shell
pixi run check-spelling
```

Or, from an environment with `codespell` installed, simply

```shell
codespell
```

To fix spelling errors in all files, use `-w` (`--write-changes`). When run in this way, the tool will modify the file in place. If unresolvable errors are encountered, these are written to standard output and must be manually fixed before attempting to rerun the tool.

**Note**: Spell checking by codespell may make unwanted changes (for example, a variable name in source code). As a result, you should check the `codespell` changes. codespell can be forced to leave a particular word unchanged by adding it to the `.codespell.ignore` file.

### Fortran formatting

Fortran source code can be formatted with [fprettify](https://github.com/pseewald/fprettify), specifying the [MODFLOW 6 fprettify configuration](.fprettify.yaml). The `fprettify` package is included in the Conda `environment.yml` and the Pixi `pixi.toml` files and can be run directly, via Pixi, or via [VSCode](.vscode/README.md) tasks.

For instance, to format a single file:

```shell
fprettify -c .fprettify.yaml ./utils/zonebudget/src/zbud6.f90
```

When run in this way, the tool will modify the file in place and generate no output if successful. If unresolvable formatting errors are encountered (e.g. for excess line length), these are written to standard output and must be manually fixed before attempting to rerun the tool.

To check whether the repository's source files satisfy formatting requirements without making any changes:

```shell
python .github/common/check_format.py
```

or using pixi:

```shell
pixi run check-format
```

To format all files, add the `--write-changes` flag to the end of the python or pixi commands. These commands will exclude the proper files from formatting, including vendored library sources in [`src/Utilities/Libraries`](src/Utilities/Libraries/).

**Note**: as `fprettify` may shift code in unexpected ways, it is a good idea to visually check source files afterwards.


### Python formatting

Python code and scripts can be formatted with [ruff](https://docs.astral.sh/ruff/), specifying the [MODFLOW 6 ruff configuration](.github/common/ruff.toml). The `ruff` package is included in the Conda `environment.yml` and Pixi `pixi.toml` files and can be run directly, via Pixi, or via [VSCode](.vscode/README.md) tasks.

For instance, to format a single file:

```shell
ruff format autotest/test_gwe_cnd.py
```

When run in this way, `ruff` will modify the file in place and generate no output if successful. If unresolvable formatting errors are encountered, these are written to standard output and must be manually fixed before attempting to rerun the tool.

To check whether the repository's python code and scripts satisfy formatting requirements without making any changes:

```shell
ruff format --check .
```

or using pixi:

```shell
pixi run check-python-format
```

To format all files, remove the `--check` flag from the python command or run the pixi command:

```shell
pixi run fix-python-format
```

### Python linting

Linting is the automated checking of source code for programmatic and stylistic errors. python code and scripts can be linted with [ruff](https://docs.astral.sh/ruff/), specifying the [MODFLOW 6 ruff configuration](.github/common/ruff.toml). The `ruff` package is included in the Conda `environment.yml` and Pixi `pixi.toml` files and can be run directly, via Pixi, or via [VSCode](.vscode/README.md) tasks.

For instance, to lint a single file:

```shell
ruff check --fix autotest/test_gwe_cnd.py
```

When run in this way, `ruff` will modify the file in place and generate no output if successful. If unresolvable formatting errors are encountered, these are written to standard output and must be manually fixed before attempting to rerun the tool.

To check whether the repository's python code and scripts satisfy linting requirements without making any changes:

```shell
ruff check .
```

or using pixi:

```shell
pixi run check-python-lint
```

To format all files, add the `--fix` flag to the python command or pixi command. Alternatively with pixi run:

```shell
pixi run fix-python-lint
```

## Testing

MODFLOW 6 unit tests are written in Fortran with [`test-drive`](https://github.com/fortran-lang/test-drive).

MODFLOW 6 integration tests are written in Python with [`pytest`](https://docs.pytest.org/en/7.1.x/). Integration testing dependencies are included in Pixi and Conda environments.

**Note:** the entire test suite should pass before a pull request is submitted. Tests run in GitHub Actions CI and a PR can only be merged with passing tests. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for more information.

### Configuring a test environment

Before running tests, there are a few steps to complete. Most importantly, the local development version of MODFLOW 6 must be built, e.g. with Meson as described above.

The `autotest/build_exes.py` script is provided as a shortcut to rebuild local binaries. It can be invoked as a standard Python script or with Pytest. By default, binaries are placed in the `bin` directory relative to the project root, as in the Meson commands described above. To change the location of the binaries, use the `--path` option.

#### Configuring unit tests

Unit tests are [driven with Meson](https://mesonbuild.com/Unit-tests.html). A small number of Meson-native tests are defined in the top-level `meson.build` file to check that MODFLOW 6 has installed successfully. These require no additional configuration.

Additional Fortran unit tests are defined with [`test-drive`](https://github.com/fortran-lang/test-drive) in the `autotest/` folder, with test files named `Test*.f90`. If Meson fails to find the `test-drive` library via `pkg-config`, these will be skipped.

To install `test-drive`:

1. Clone the `test-drive` repository
2. Setup/build with Meson, e.g. in a Unix shell from the `test-drive` project root:

```shell
meson setup builddir --prefix=$PWD --libdir=lib
meson install -C builddir
```

3. Add `<test-drive project root>/lib/pkgconfig` to the `PKG_CONFIG_PATH` environment variable.
4. To confirm that `test-drive` is detected by `pkg-config`, run `pkg-config --libs test-drive`.

Meson should now detect the `test-drive` library when building MODFLOW 6.

**Note:** the `test-drive` source code is not yet compatible with recent versions of Intel Fortran, building with `gfortran` is recommended.

See the [Running unit tests](#running-unit-tests) section for instructions on running unit tests.

#### Configuring integration tests

A few more tasks must be completed before integration testing:

- install MODFLOW-related executables
- ensure FloPy packages are up to date
- install MODFLOW 6 example/test models

As mentioned above, binaries live in the `bin` subdirectory of the project root. This directory is organized as follows:

- local development binaries in the top-level `bin`
- binaries rebuilt in development mode from the latest MODFLOW 6 release in `bin/rebuilt/`
- related programs installed from the [executables distribution](https://github.com/MODFLOW-ORG/executables/releases) in `bin/downloaded/`

##### Rebuilding release binaries

Tests require the latest official MODFLOW 6 release to be compiled in develop mode with the same Fortran compiler as the development version. A number of binaries distributed from the [executables repo](https://github.com/MODFLOW-ORG/executables) must also be installed. The script `autotest/get_exes.py` does both of these things. It can be run from the project root with:

```shell
pixi run get-exes
```

Alternatively, from the `autotest/` directory:

```shell
pytest get_exes.py
```

As above, binaries are placed in the `bin` subdirectory of the project root, with nested `bin/downloaded` and `bin/rebuilt` subdirectories containing the rebuilt latest release and downloaded binaries, respectively.

##### Updating FloPy packages

FloPy packages should be regenerated from DFN files before running tests for the first time or after definition files change. This can be done with the `autotest/update_flopy.py` script, which wipes and regenerates package classes for the FloPy installed in the Python environment.

**Note:** if you've installed an editable local version of FloPy from source, running this script can overwrite files in your repository.

There is a single optional argument, the path to the folder containing definition files. By default DFN files are assumed to live in `doc/mf6io/mf6ivar/dfn`, making the following functionally identical:

```shell
pixi run update-flopy
```

which uses the default dfn path. Or the location of the definition files can be explitily defined using:

```shell
pixi run update-flopy doc/mf6io/mf6ivar/dfn
```

Alternatively, run `python update_flopy.py` directly from `autotest/`.

##### Updating Fortran definitions

Any time a MODFLOW 6 input definition file (dfn) has been changed internal MODFLOW 6 Fortran definitions should be updated as well. This can be accomplished locally by running `utils/idmloader/scripts/dfn2f90.py` and then recompiling. This script will update the appropriate input package Fortran definition files if the dfn change is relevant to input processing. Updated Fortran definition files should accompany related dfn file changes when creating a pull request.

```shell
cd utils/idmloader/scripts
python dfn2f90.py
```

or using pixi:

```shell
pixi run update-fortran-definitions
```


##### Installing external models

Some autotests load models from external repositories:

- [`MODFLOW-ORG/modflow6-testmodels`](https://github.com/MODFLOW-ORG/modflow6-testmodels)
- [`MODFLOW-ORG/modflow6-largetestmodels`](https://github.com/MODFLOW-ORG/modflow6-largetestmodels)
- [`MODFLOW-ORG/modflow6-examples`](https://github.com/MODFLOW-ORG/modflow6-examples)

By default, the test framework will test MODFLOW 6 against these models as accessed via the [MODFLOW devtools models API](https://modflow-devtools.readthedocs.io/en/latest/md/models.html). It may be necessary to test MODFLOW 6 against models on the local filesystem. See the [MODFLOW devtools documentation](https://modflow-devtools.readthedocs.io/en/latest/md/install.html#installing-external-model-repositories) for instructions to clone and install external model repositories.

### Running tests

MODFLOW 6 has two kinds of tests: Fortran unit tests, driven with Meson, and Python integration tests, driven with Pytest.

#### Running unit tests

Unit tests must be run from the project root. To run unit tests in verbose mode:

```shell
meson test -C builddir
```
or using pixi:

```shell
pixi run test builddir
```

Unit tests can be selected by module name (as listed in `autotest/tester.f90`). For instance, to test the `ArrayHandlersModule`:

```shell
meson test -C builddir --verbose ArrayHandlers
```

or using pixi:

```shell
pixi run test builddir --verbose ArrayHandlers
```

To run a test module in the `gdb` debugger, just add the `--gdb` flag to the test command.

#### Running integration tests

Integration tests must be run from the `autotest/` folder if invoked with `pytest` directly &mdash; the Pixi `autotest` task can be invoked from the project root.

To run tests in parallel:

```shell
cd autotest/
pytest -v -n auto  # from autotest/
```

or using pixi:

```shell
pixi run autotest
```

The Pixi `autotest` task includes options to run tests in parallel, show test runtimes, and save failed test results in `autotest/.failed/`.

**Note:** The `-n` option accepts an integer argument for the number of parallel processes. If the value `auto` is provided, `pytest-xdist` will use one worker per available processor.

##### Selecting tests with markers

Markers can be used to select subsets of tests. Markers provided in `pytest.ini` include:

- `slow`: tests that take longer than a few seconds to complete
- `external`: tests that use models in external repositories
- `large`: tests that use large models
- `regression`: tests comparing results from multiple versions

Markers can be used with the `-m <marker>` option, and can be applied in boolean combinations with `and`, `or` and `not`. For instance, to run fast tests in parallel, excluding regression tests:

```shell
pytest -v -n auto -m "not slow and not regression"
```

The `--smoke` (short `-S`) flag, provided by `modflow-devtools` is an alias for the above:

```shell
pytest -v -n auto -S
```

or using pixi:

```shell
pixi run autotest -S
```

[Smoke testing](https://modflow-devtools.readthedocs.io/en/latest/md/markers.html#smoke-testing) is a form of integration testing which aims to test a decent fraction of the codebase quickly enough to run often during development.

Tests using models from external repositories can be selected with the `external` marker:

```shell
pixi run autotest -m "external"
```

By default, these will run against test models pulled from the GitHub repositories. To run the tests against local models, use `--models-path` once or more to specify directories to search for model input files. For instance, to test MF6 models from the test models repository:

```shell
pixi run autotest -m "external" --models-path /path/to/modflow6-testmodels/mf6
```

By defaultk, only MF6 models are found. To test the mf5to6 converter with mf2005 models in the same repository, relax the namefile search pattern:

```shell
pixi run autotest -m "external" --models-path /path/to/modflow6-testmodels/mf5to6 --namefile-pattern "*.nam"
```

Large test models are excluded from commit-triggered CI and only run on GitHub Actions nightly. To run the models from a local clone of the repository:

```shell
pixi run autotest -m "external" --models-path /path/to/modflow6-largetestmodels
```

Tests load external models from fixtures provided by `modflow-devtools`. External model tests can be selected by model or simulation name, or by packages used. See the [`modflow-devtools` documentation](https://modflow-devtools.readthedocs.io/en/latest/md/fixtures.html#filtering) for usage examples. Note that filtering options only apply to tests using external models, and will not filter tests defining models in code &mdash; for that, the `pytest` built-in `-k` option may be used.

### Writing tests

#### Writing unit tests

To add a new unit test:

- Add a file containing a test module, e.g. `TestArithmetic.f90`, to the `autotest/` folder.

```fortran
module TestArithmetic
  use testdrive, only : error_type, unittest_type, new_unittest, check, test_failed
  implicit none
  private
  public :: collect_arithmetic
contains
  
  subroutine collect_arithmetic(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [new_unittest("add", test_add)]
  end subroutine collect_arithmetic

  subroutine test_add(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, 1 + 1 == 2, "Math works")
    if (allocated(error)) then
      call test_failed(error, "Math is broken")
      return
    end if
  end subroutine test_add
end module TestArithmetic
```

- Add the module name to the list of `tests` in `autotest/meson.build`, omitting the leading "Test".

```fortran
tests = [
  'Arithmetic',
]
```

- Add a `use` statement for the test module in `autotest/tester.f90`, and add it to the array of `testsuites`.

```fortran
use TestArithmetic, only: collect_arithmetic
...
testsuites = [ &
  new_testsuite("Arithmetic", collect_arithmetic), &
  new_testsuite("something_else", collect_something_else) &
]
```

- Rebuild with Meson from the project root, e.g. `meson install -C builddir`. The test should now be picked up when `meson test...` is next invoked.

#### Writing integration tests

Integration tests should ideally follow a few conventions for easier maintenance:

- Use temporary directory fixtures. Tests which write to disk should use `pytest`'s built-in `tmp_path` fixtures or one of the [keepable temporary directory fixtures from `modflow-devtools`](https://modflow-devtools.readthedocs.io/en/latest/md/fixtures.html#keepable-temporary-directories). This prevents tests from polluting one another's state.

- Use markers for convenient (de-)selection:
  - `@pytest.mark.slow` if the test doesn't complete in a few seconds (this preserves the ability to quickly [`--smoke` test](https://modflow-devtools.readthedocs.io/en/latest/md/markers.html#smoke-testing)
  - `@pytest.mark.external` if the test relies on external model repositories
  - `@pytest.mark.regression` if the test compares results from different versions

**Note:** If all three external model repositories are not installed as described above, some tests will be skipped. The full test suite includes >750 cases. All must pass before changes can be merged into this repository.

##### Test framework

A framework has been developed to streamline common testing patterns. The [`TestFramework`](autotest/framework.py) class, defined in `autotest/framework.py`, is used by most test scripts to configure, run and evaluate one or more MF6 simulations, optionally in comparison with another simulation or model.

Generally, the recommended pattern for a test script is:

```python
import ...

cases = ["a", "b", ...]
variable = [1., 0., ...]
expected = [-1., -1.1, ...]

def build_models(idx, test):
  v = variable[idx]
  ...

def check_output(idx, test):
  e = expected[idx]
  ...

def plot_output(idx, test):
  import matplotlib.pyplot as plt
  ...

@pytest.mark.parametrize("idx, name", enumerate(cases))
def test_mf6model(idx, name, function_tmpdir, targets, plot):
    test = TestFramework(
        name=name,
        workspace=function_tmpdir,
        targets=targets,
        build=lambda t: build_models(idx, t),
        check=lambda t: check_output(idx, t),
        plot=lambda t: plot_output(idx, t) if plot else None,
        compare=None,
    )
    test.run()
```

The framework has three hooks:

- `build`: construct one or more MF6 simulations and/or non-MF6 models with FloPy
- `check`: evaluate simulation/model output
- `plot`: make one or more plots of simulation/model output

A test script conventionally contains one or more test cases, fed to the test function as `idx, name` pairs. `idx` can be used to index parameter values or expected results for a specific test case. The test case `name` is useful for model/subdirectory naming, etc.

The framework will not run an unknown program. The path to any program under test (or used for a comparison) must be registered in the `targets` dictionary. Keys are strings. See `autotest/conftest.py` for the contents of `targets` &mdash; naming follows the [executables distribution](https://github.com/MODFLOW-ORG/executables).

The `.run()` function

1. builds simulations/models
2. runs simulations/models
3. compares simulation/model outputs
4. checks outputs against expectations

A `compare` parameter may be provided on initialization, which enables comparison of outputs against another program or the latest official release of MF6. The following values are supported:

- `None`: disables comparison &mdash; the test simply runs/evaluates any registered simulations/models without comparing results
- `auto`: attempt to detect the comparison type from contents of test workspace, otherwise skipping comparison
- `mf6_regression`: compare results against the latest official release rebuilt in develop mode
- `mf6`, `mf2005`, `mfnwt`, or `mflgr`: compare with results from the selected program &mdash; a corresponding model must be provided in `build_models()`

After running the reference and comparison models, the framework will try to find correspondingly named output files to compare &mdash; comparison logic may need adjustment when writing tests for new packages or models.

## Generating makefiles

Run `build_makefiles.py` in the `distribution/` directory after adding, removing, or renaming source files. This script uses [Pymake](https://github.com/modflowpy/pymake) to regenerate makefiles. For instance:

```shell
cd distribution/
python build_makefiles.py
```

or using pixi:

```shell
pixi run build-makefiles
```

### Updating extra and excluded files

If the utilities located in the `utils` directory (e.g., `mf5to6` and `zbud6`) are affected by changes to the modflow6 `src/` directory (such as new or refactored source files), then the new module source file should also be added to the utility's `utils/<util>/pymake/extrafiles.txt` file. This file informs Pymake of source files living outside the main source directory, so they can be included in generated makefiles.

Module dependencies for features still under development should be added to `excludefiles.txt`. Source files listed in this file will be excluded from makefiles generated by Pymake. Makefiles should only include the source files needed to the build officially released/supported features.

### Testing makefiles

Makefile generation and usage can be tested from the `distribution` directory by running the `build_makefiles.py` script with Pytest:

```shell
pytest -v build_makefiles.py
```

**Note**: `make` is required to test compiling MODFLOW 6 with makefiles. If `make` is not discovered on the system path, compile tests will be skipped.

Makefiles may also be tested manually by changing to the appropriate `make` subdirectory (of the project root for MODFLOW 6, or inside the corresponding `utils` subdirectory for the zonebudget or converter utilities) and invoking `make` (`make clean` may first be necessary to remove previously created object files).

### Installing `make` on Windows

On Windows, it is recommended to generate and test makefiles from a Unix-like shell rather than PowerShell or Command Prompt. Make can be installed via [Conda](https://anaconda.org/conda-forge/make) or [Chocolatey](https://community.chocolatey.org/packages/make). Alternatively, it is included with [mingw](https://sourceforge.net/projects/mingw/), which is also available from [Chocolatey](https://community.chocolatey.org/packages/mingw).

#### Using Conda from Git Bash

To use Conda from Git Bash on Windows, first run the `conda.sh` script located in your Conda installation's `/etc/profile.d` subdirectory. For instance, with Anaconda3:

```shell
. /c/Anaconda3/etc/profile.d/conda.sh
```

Or Miniconda3:

```shell
. /c/ProgramData/miniconda3/etc/profile.d/conda.sh
```

After this, `conda` commands should be available.

This command may be added to a `.bashrc` or `.bash_profile` file in your home directory to permanently configure Git Bash for Conda.

## Branching model

This section documents MODFLOW 6 branching strategy and other VCS-related procedures.

### Overview

This project follows the [git flow](https://nvie.com/posts/a-successful-git-branching-model/): development occurs on the `develop` branch, while `master` is reserved for the state of the latest release. Development PRs are typically squashed to `develop` to avoid merge commits. At release time, release branches are merged to `master`, and then `master` is merged back into `develop`.

### Managing long-lived branches

When a feature branch takes a long time to develop, it is easy to become out of sync with the develop branch.  Depending on the situation, it may be advisable to periodically squash the commits on the feature branch and rebase the change set with develop.  The following approach for updating a long-lived feature branch has proven robust.

In the example below, the feature branch is assumed to be called `feat-xyz`.

#### Backup

Begin by creating a backup copy of the feature branch in case anything goes terribly wrong.

```
git checkout feat-xyz
git checkout -b feat-xyz-backup
git checkout feat-xyz
```

#### Squash

Next, consider squashing commits on the feature branch.  If there are many commits, it is beneficial to squash them before trying to rebase with develop.  There is a nice article on [squashing commits into one using git](https://www.internalpointers.com/post/squash-commits-into-one-git), which has been very useful for consolidating commits on a long-lived modflow6 feature branch.

A quick and dirty way to squash without interactive rebase (as an alternative to the approach described in the article mentioned in the preceding paragraph) is a soft reset followed by an amended commit. First making a backup of the feature branch is strongly recommended before using this approach, as accidentally typing `--hard` instead of `--soft` will wipe out all your work.

```
git reset --soft <first new commit on the feature branch>
git commit --amend -m "consolidated commit message"
```

Once the commits on the feature branch have been consolidated, a force push to origin is recommended.  This is not strictly required, but it can serve as an intermediate backup/checkpoint so the squashed branch state can be retrieved if rebasing fails.  The following command will push `feat-xyz` to origin.

```
git push origin feat-xyz --force
```

The `--force` flag's short form is `-f`.

#### Rebase

Now that the commits on `feat-xyz` have been consolidated, it is time to rebase with develop.  If there are multiple commits in `feat-xyz` that make changes, undo them, rename files, and/or move things around in subsequent commits, then there may be multiple sets of merge conflicts that will need to be resolved as the rebase works its way through the commit change sets.  This is why it is beneficial to squash the feature commits before rebasing with develop.

To rebase with develop, make sure the feature branch is checked out and then type:

```
git rebase develop
```

If anything goes wrong during a rebase, there is the `rebase --abort` command to unwind it.

If there are merge conflicts, they will need to be resolved before going forward.  Once any conflicts are resolved, it may be worthwhile to rebuild the MODFLOW 6 program and run the smoke tests to ensure nothing is broken.  

At this point, you will want to force push the updated feature branch to origin using the same force push command as before.

```
git push origin feat-xyz --force
```

#### Cleanup

Lastly, if you are satisfied with the results and confident the procedure went well, then you can delete the backup that you created at the start.

```
git branch -d feat-xyz-backup
```

This process can be repeated periodically to stay in sync with the develop branch and keep a clean commit history.

## Deprecation policy

To deprecate a MODFLOW 6 input/output option in a DFN file:

- Add a new `deprecated x.y.z` attribute to the appropriate variable in the package DFN file, where `x.y.z` is the version the deprecation is introduced. Mention the deprecation prominently in the release notes.
- If support for the deprecated option is removed (typically after at least 2 minor or major releases or 1 year), add a new `removed x.y.z` attribute to the variable in the DFN file, where `x.y.z` is the version in which support for the option was removed. The line containing `deprecated x.y.z` should not be deleted. Mention the removal prominently in the release notes.
- Deprecated/removed attributes are not removed from DFN files but remain in perpetuity. The `doc/mf6io/mf6ivar/deprecations.py` script generates a markdown deprecation table which is converted to LaTeX by `doc/ReleaseNotes/mk_deprecations.py` for inclusion in the MODFLOW 6 release notes. Deprecations and removals should still be mentioned separately in the release notes, however.

### Finding deprecations

To search for deprecations and removals in DFN files on a system with `git` and standard Unix commands available:

```shell
git grep 'deprecated' -- '*.dfn' | awk '/^*.dfn:deprecated/'
```

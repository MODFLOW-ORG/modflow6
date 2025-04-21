# Developing extended MODFLOW 6

This document describes how to set up your build environment for developing and testing the extended version of MODFLOW. It further builds on the instructions given in [DEVELOPER.md](DEVELOPER.md) so make sure to read those first.

---
**DISCLAIMER**

*Expectations on platform compatibility*

The serial version of the MODFLOW 6 program has had no external dependencies and is traditionally available for a variety of platforms (Windows, GNU/linux, macOS) and compatible with the mainstream Fortran compilers (gfortran, ifort). The extended version comes with dependencies on third party components, most notably the MPI, PETSc and NetCDF libraries. While the goal is a continued support of the above mentioned configurations, this has become more challenging and can generally not be guaranteed. To assist developers as well as end users who are planning to compile the code themselves, a list of successfully tested build configurations will be included in this document.

---

<!-- TODO: insert TOC tags here -->

## Introduction

The design philosophy has been to maintain MODFLOW as a single codebase and have it compile to either a serial or extended program. The former continues to be a lightweight, highly compatible code which does not require external dependencies other than those provided by the standard compiler libraries. The latter has three distinct dependencies on 3rd party libraries: MPI, PETSc and NetCDF, as described below. Note that extended capabilities are a true extension and the executable will in all cases be capable of *serial* execution with equivalent results.

## Prerequisites

The extended version of MODFLOW 6 requires the Message Passing Interface (MPI), the Portable, Extensible Toolkit for Scientific Computation (PETSc - pronounced PET-see (/ˈpɛt-siː/)) and the NetCDF for array-oriented scientific data libraries.

### MPI

The extended version of MODFLOW 6 uses MPI to synchronize data between processes when configured to run in parallel mode. There are a couple of implementations of the MPI standard available. Their applicability usually depends on the platform that is used:

- Open MPI: https://www.open-mpi.org/
- MPICH: https://www.mpich.org/
- Intel MPI: https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html
- Microsoft MPI: https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi

On Linux and macOS, if you haven't installed a MPI framework yet, your best bet is to automatically download the MPI implementation if building PETSc from source using the option flag `--download-openmpi` or `--download-mpich` and the `configure` approach described below. Alternatively, you can install OpenMPI or MPICH using a package manager for your OS. Tested configurations are listed in the table below.

In addition to compiling, the MPI toolset is also required to run a parallel simulation. The implementations above all come with `mpiexec` (or `mpiexec.exe`) to start an executable in parallel.

### PETSc

The PETSc library is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations:

https://petsc.org/release/

The PETSc library (version 3.16 or higher) is used by MODFLOW for its parallel linear solver capabilities and the distributed data formats (vectors and matrices) that go along with it. Parallel PETSc uses MPI internally as well, so setting up this library should typically be coordinated with the installation of the MPI library. A lot of obscure things can happen if the binaries are not compatible, so in general it is a good strategy to compile MPI, PETSc, and MODFLOW with the same compiler toolchain or install MPI and PETSc using a package manager and built with the same compiler being used to compile MODFLOW.

## Compiling MPI and PETSC from source

The PETSc website gives details on a large number of configurations, depending on the target platform/OS, and many different ways to configure/make/install the library: https://petsc.org/release/install/. Building on Windows is notoriously challenging and discouraged by the PETSc development team. On Linux, however, PETSc can be installed (configure/make/install) by executing the following command

```
$ ./configure --download-openmpi --download-fblaslapack
$ make all
```

in a terminal open in the root directory of your PETSc download

## Compiling NetCDF from source

See, for example, https://docs.unidata.ucar.edu/netcdf-c/4.9.2/building\_netcdf\_fortran.html and https://docs.unidata.ucar.edu/nug/current/getting\_and\_building\_netcdf.html#building for information related to building NetCDF-Fortran and prerequisites in a Unix/Linux environment.

## Using a package manager to install MPI and PETSc

Use of a package manager can simplify the process of building the extended version of MODFLOW 6.

### MacOS

[OpenMPI](https://formulae.brew.sh/formula/open-mpi) and [PETSc](https://formulae.brew.sh/formula/petsc) are available on Homebrew for Intel and Apple Silicon (M1). Both of these depend on [gcc 13.1.0](https://formulae.brew.sh/formula/gcc). [NetCDF Fortran](https://formulae.brew.sh/formula/netcdf-fortran) and related dependencies are also available on Homebrew for Intel and Apple Silicon (M1). [pkg-config](https://formulae.brew.sh/formula/pkg-config) should also be installed from Homebrew, if not already installed, so that Meson will be able to resolve the installation location of MPI and PETSc.

### Ubuntu

OpenMPI, PETSc and NetCDF Fortran are available for a variety of Ubuntu versions using the Advanced Packaging Tool (apt).

### Windows

Under evaluation.

## Using pkg-config to check your PETSc installation

Eventually, the MODFLOW build process has to resolve the installation location of all external dependencies. The pkg-config tool (https://en.wikipedia.org/wiki/Pkg-config) can be used to take care of that.

```
pkg-config --libs petsc
```

If PETSc was build from source, you can check the contents of the folder

```
$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig/
```

and confirm that there are one or more `*.pc` files in there. A similar `pkgconfig/` folder has to present for the MPI installation that was used. For example, for Open MPI on WSL2 this folder is `/lib/x86_64-linux-gnu/pkgconfig/`.

To connect everything, both of these folder paths have to be added to the `PKG_CONFIG_PATH` variable so that the `pkg-config` executable can resolve the installed libraries.

## Building the extended version of MODFLOW 6

The primary build system for MODFLOW is Meson (https://mesonbuild.com/). The `meson.build` script takes an additional argument to activate an extended build of the software. E.g for building and installing an extended release version:

```
meson setup builddir -Ddebug=false -Dextended=true \
                     --prefix=$(pwd) --libdir=bin
meson install -C builddir
meson test --verbose --no-rebuild -C builddir
```

Note that changing the option flags in the `meson setup` command requires the flag `--reconfigure` to reconfigure the build directory. If the `PKG_CONFIG_PATH` was set as described above, the linking to PETSc and MPI is done automatically.

It's always a good idea to check your extended MODFLOW executable to confirm that it is successfully linked against the external dependencies. You can use the command line tools `ldd` (Linux), `otool` (macOS), or `Dependencies.exe` (Windows, https://github.com/lucasg/Dependencies) to do that. In the list of dependencies, you should be able to identify `libpetsc`, `libmpi`, `libnetcdff`, `libnetcdf`, and `libhdf5` for extended builds.

The other build systems in the MODFLOW project (MS Visual Studio, `pymake`, `Makefile`) continue to be supported for *serial* builds only. `pymake` uses the `excludefiles.txt` to ignore those files that can only be build when MPI and PETSc are present on the system. In MS Visual Studio these same files are included in the solution but not in the build process.

---

**IMPORTANT**

*Don't use MPI and PETSc directly in your code*

Extended MODFLOW was designed to have all third party functionality (MPI, PETSc and the NetCDF Fortran API) made available through the framework. Developers of models and packages **should not** directly call these libraries and change the set of excluded files described above. If you feel you need to include MPI or PETSc functionality in your code (e.g. you want to `use mpi` in your source file), contact the MODFLOW development team on how to best proceed.

---

## Testing the extended version of MODFLOW 6

Extended MODFLOW can be tested using the same test framework as the serial program, with just a few modifications. To run a test inside the `autotest` folder in parallel mode, make sure to add a marker `@pytest.mark.parallel` so that the test is only executed in the Continuous Integration when running a configuration with an extended build of MODFLOW. Similarly, adding the marker `@pytest.mark.netcdf` to a test ensures that the test will only be executed when intended by specifying NetCDF test mode.

The `TestSimulation` object that is being run from the framework should be configured for parallel run mode with the flag `parallel=True` and the number of processes `ncpus=...`. As an example, see the test case in `autotest/test_par_gwf01.py`, which can be run with

```
$ pytest -s --parallel test_par_gwf01.py
```

To run a test in NetCDF test mode, for example `autotest/test_netcdf_gwe_cnd.py`:

```
$ pytest --netcdf test_netcdf_gwe_cnd.py
```

Running without the `--parallel` and `--netcdf` flags will simply skip these tests.

## Debugging

The most straightforward way to debug a parallel simulation is to start a run and have it pause to attach the debugger(s). Make sure that the MODFLOW executable was compiled with `-Ddebug=true`. In parallel mode the program uses the PETSc solver and a configuration file `.petscrc` should be present in the same folder as the simulation's `mfsim.nam`. In that PETSc resource file, you should add the following option:

```
-wait_dbg
```

telling MODFLOW to pause immediately after startup. This will give you time to attach one or multiple debuggers to the processes. Then start the parallel program, for example on two cores:

```
mpiexec -np 2 mf6 -p
```

In the process explorer you should now see 2 processes called `mf6` or `mf6.exe`. On the prompt where the command was executed, MODFLOW waits for input:

```
$ Hit enter to continue...
```

In order to truly debug in parallel, i.e. step through the instructions side-by-side, you will need to start two instances of the debugger and attach them. The following section describes how to do that with VSCode.

### Debugging with VSCode

In VSCode parallel debugging is easiest done by duplicating the development environment. First, make sure that you have set up your environment to build an extended version in debug mode. The VSCode launch.json should contain an entry to attach to a running process:

```
  {
    "name": "Attach to ...",
    "type": "cppdbg",
    "request": "attach",
    "processId": "${command:pickProcess}",
    "program": "<path_to>/mf6",
    "MIMode": "gdb",
    "miDebuggerPath": "<path_to>/gdb",
    "setupCommands": [
      {
        "description": "Enable pretty-printing for gdb",
        "text": "-enable-pretty-printing",
        "ignoreFailures": true
      }
    ]
  }
```

After building extended MODFLOW, press `Ctrl+Shift+p` to execute *Workspaces: Duplicate As Workspace in New Window*. This will open a second VSCode window, identical to the first. Starting the debug process and selecting *"Attach to ..."* opens a process selection window with the processes started from the `mpiexec` command described above. Select both, each from their own instance of the VSCode program. Now you can put breakpoints in the code, "Hit enter to continue" on the command prompt, and step through the parallel processes side-by-side.

---
**TIP**

Make sure that you work with gdb versions >= 10. We have found that earlier versions are only partially compatible with the VSCode debugging interface and crash when inspecting the data of Fortran derived types

---

## Compatibility

Parallel MODFLOW (a pre-extended build without NetCDF) has been built successfully with the following configurations:

| Operating System                    | Toolchain                 | MPI               | PETSc               | Package Manager |
|-------------------------------------|---------------------------|-------------------|---------------------|-----------------|
| MS Windows                          | ?                         | ?                 | ?                   | NA              |
| WSL2 (Ubuntu 20.04.5)               | gcc 9.4.0                 | OpenMPI 4.0.3     | 3.18.2              | NA              |
| macOS 12.6.3                        | gcc 9.5.0                 | OpenMPI 4.1.4     | 3.18.5              | NA              |
| macOS 12.6.6                        | gcc 13.1.0                | OpenMPI 4.1.5     | 3.19.1              | Homebrew        |
| Ubuntu 22.04                        | gcc 9.5.0                 | OpenMPI 4.1.4     | 3.18.5              | NA              |
| Ubuntu 22.04 ARM64                  | gcc 11.4.0                | OpenMPI 4.1.5     | 3.19.3              | apt             |
| Ubuntu 22.04 ARM64                  | gcc 9.5.0                 | MPICH 3.4.1       | 3.15.5              | NA              |
| Ubuntu 22.04 ARM64                  | gcc 12.3.0                | MPICH 4.1.1       | 3.19.6              | NA              |
| Ubuntu 22.04 ARM64                  | gcc 12.3.0                | MPICH 4.1.1       | 3.20.0              | NA              |
| SUSE Linux Enterprise Server 15 SP2 | intel 19.1.0.166 20191121 | CRAY-MPICH 7.7.19 | CRAY-PETSC 3.14.5.0 | NA              |
| Red Hat Enterprise Linux 8.7        | intel 2021.10.0 20230609  | CRAY-MPICH 8.1.26 | 3.15.5              | NA              |

The most up-to-date configurations are available in the GitHub CI script: `.github/workflows/ci.yml` under the task `extended_test`. These are being tested upon every change to the `develop` branch of MODFLOW.

To improve support, we kindly ask you to share your experience with building and running extended MODFLOW and report back if you have a successful setup that is not in this table.

## Known issues

### Building PETSc on Ubuntu 22.04 with MPICH and GNU compilers

Versions of PETSc that use MPICH 3.4 (v3.14, v3.15, v3.16)) must be built with gcc-9 or earlier. Versions of PETSc that use MPICH 4.1 (v3.17 or newer) can be built with newer versions of the gcc compiler (gcc-11, gcc-12, etc.).

Meson does not correctly load the Fortran compiler flags from the `mpich.pc` package configuration file in the `$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig` directory. To overcome this issue, make a copy of `mpich.pc` and name it `mpichfort.pc`. Then determine the appropriate Fortran flags using

```
$PETSC_DIR/$PETSC_ARCH/bin/mpifort -show
```

which will return something like

```
$ linux-real-gcc12.3.0-3.20.0/bin/mpifort -show
gfortran -fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O2 -fallow-argument-mismatch -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -L/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -lmpifort -Wl,-rpath -Wl,/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -Wl,--enable-new-dtags -lmpi
```

Copy the returned Fortran flags and replace the `Libs:` and `Cflags:` attributes in  the `mpichfort.pc` file. Also modify the `Name:` attribute to `mpichfort`. The modified `mpichfort.pc` file should look something like

```
# this gives access to the mpich header files
prefix=/media/psf/Development/petsc/linux-real-gcc9.5.0-3.20.0
exec_prefix=${prefix}
libdir=/media/psf/Development/petsc/linux-real-gcc9.5.0-3.20.0/lib
includedir=${prefix}/include

Name: mpichfort
Description: High Performance and portable MPI
Version: 4.1.2
URL: http://www.mcs.anl.gov/research/projects/mpich
Requires:
Libs: -fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O2 -fallow-argument-mismatch -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -L/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -lmpifort -Wl,-rpath -Wl,/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -Wl,--enable-new-dtags -lmpi 
Cflags:  -fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O2 -fallow-argument-mismatch -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -I/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include -L/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -lmpifort -Wl,-rpath -Wl,/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib -Wl,--enable-new-dtags -lmpi

# pkg-config does not understand Cxxflags, etc. So we allow users to
# query them using the --variable option

cxxflags= -Wno-lto-type-mismatch -Wno-psabi -O2 -std=gnu++17 -fPIC -I${includedir}
fflags=-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O2 -I${includedir}
fcflags=-fPIC -ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -O2 -I${includedir}
```

The `/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/include` and `/media/psf/Development/petsc/linux-real-gcc12.3.0-3.20.0/lib` entries in the `Libs:` and `Cflags:` attributes can be replaced with `${includedir}` and `${libdir}`, respectively, to simplify `mpichfort.pc`.

# MODFLOW 6 Input Data Model
MODFLOW 6 reads simulation input from text and binary input files.  Traditionally, each model package has been responsible for reading its own input.  Although there are low-level utilities for reading arrays, and parsing input blocks, these utilities are called from dedicated read routines for each package.  To support new types of input, such as from NetCDF files or through the Application Programming Interface, a new effort is underway to implement a comprehensive Input Data Model (IDM) for MODFLOW 6.  Implementation of an Output Data Model may follow this effort.

- [Overview](#overview)
- [Framework](#framework)
- [Design](#design)
- [Integration](#integration)

## Overview
The MODFLOW 6 IDP (Input Data Processor) is a subsystem meant to generically read and store simulation input data from more than one type of input source.  It is built upon existing MODFLOW 6 capabilities, specifically input parameter descriptions currently stored in \*.dfn (definition) files, and the Memory Manager.  The parameter descriptions provide IDP with input block and parameter attributes needed to create memory and store input data.  The Memory Manager provides a globally accessible way to create, update, access and remove input (and other) data.  Existing components that read a traditional ASCII input file can be updated to source their input from the Memory Manager.  This type of update is a significant step towards supporting sources of input data other than the traditional ASCII files. 

## Framework
MODFLOW 6 \*.dfn files are used pre-compile time to generate fortran source files containing a subset of the parameter and block attribute information necessary to process input for the simulation. A single definition file is converted into one fortran source file that defines the component parameters and blocks and organizes them into lists.  This conversion from definition file to fortran source file is currently managed by the [dfn2f90.py](scripts/dfn2f90.py) script.  This script also creates framework fortran source files with routines for generically accessing package definitions by component.  A package or other component intended to be updated and integrated with IDM must be added to the `utils/idmloader/dfns.txt` file, which the `dfn2f90.py` script consumes.  This process is described below.

## Design
Input data stored in the Memory Manager use the special identifier ```__INPUT__``` as a prefix to a memory path.  MODFLOW 6 refers to memory path prefixes such as these as a "context" and as such the collection of input data stored in the Memory Manager is referred to as the "input context".

IDM defines and implements loaders as objects to read input from some supported source and to allocate and store that input in a source independent way for the package or component.  Variations across sources are managed by the loaders.  Input retrieved by simulation packages is structured such that the source is not apparent or relevant to the package.

IDM distinguishes between static and dynamic loader objects.  Static loader objects load all pre-period input before any simulation objects are created and prepare for dynamic (stress period) data loads by creating the dynamic loader when relevant.  Static and dynamic loaders are defined per supported source by extending base types defined in [InputLoadType.f90](../../src/Utilities/Idm/InputLoadType.f90).  Top level IDM context ([IdmLoad.f90](../../src/Utilities/Idm/IdmLoad.f90)) maintains pointers to loader base types and invoke the common deferred interfaces load() (static) and rp() (dynamic) for each component.  Each source can introduce as much complexity as needed to support variations in input by introducing source specific loader types.  For example, there are currently distinct loaders for dynamic ASCII list based and array-based inputs.  These sub-loaders are allocated and managed within the context of a given source (e.g. ASCII).

## Integration
A simulation component (package, etc) can be integrated with IDM by following these steps:
- Create dfn file
- Add component to `dfns.txt` and run `dfn2f90.py`
- Update MODFLOW 6 build scripts (e.g. meson / msvs)
- Use common interfaces to source input from input context

The DFN file is already a requirement for MODFLOW 6 package integration and so is not an IDM specific requirement.  The DFN file should be added to `doc/mf6io/mf6ivar/dfn/`  and its filename added to `dfns.txt`. The generated package files are named `\*idm.f90` and located in the `src/Idm/` subdirectory by default.  To run the script:

```shell
python utils/idmloader/scripts/dfn2f90.py
```
Running this command will generate a new definition file for the package and update the IDM selector framework.  If the package also introduces a new model then a new selector framework file will also be created.  Update the MODFLOW 6 build system (meson and msvs) files with any newly generated files.

Once these files have been added to the build the package can be updated to source it's input from the input context.  In general, static input data is copied from the input context to the package (or "model") context by using the ```mem_set_value()``` interface:

```shell
call mem_set_value(this%xorigin, 'XORIGIN', this%input_mempath, found%xorigin)
```

The found instance is a convenience type generated by dfn2f90.py and should be complete and ready to use from the \*idm.f90 package module.  The interface updates the data pointer (```this%xorigin``` in the example) only if the relevant data was provided as input and will set the logical ```found%xorigin``` to TRUE when this is the case.

Dynamic data can also be copied between contexts by using the ```mem_set_value()``` interface or it can be accessed directly by setting a pointer to it.  In most cases, dynamic input memory is not reallocated during it's lifetime and such a pointer is valid for the duration of the simulation.  See CHD package code for an example (e.g. "HEAD")

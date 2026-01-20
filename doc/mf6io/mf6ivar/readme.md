# MODFLOW 6 input specification language

MODFLOW 6 accepts input via **input files** which use a custom text-based format. The MODFLOW 6 input file format is specified by **definition (DFN) files**. This document specifies the definition file format.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [DFN file overview](#dfn-file-overview)
- [Variable attributes](#variable-attributes)
  - [Reader Attribute](#reader-attribute)
- [Variable types](#variable-types)
  - [Keyword](#keyword)
  - [String](#string)
  - [Numeric](#numeric)
  - [Recarray](#recarray)
  - [Record](#record)
    - [Repeating records](#repeating-records)
  - [Keystring](#keystring)
- [Text substitutions](#text-substitutions)
- [Writing definition files](#writing-definition-files)
  - [Subpackages](#subpackages)
  - [Models](#models)
  - [Solvers](#solvers)
- [Developer tooling](#developer-tooling)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## DFN file overview

A DFN file is a text file enumerating **input variables** for a MODFLOW 6 component. An input variable is any value provided by the user to configure a MODFLOW 6 simulation. Each definition file specifies the input file expected for some MODFLOW 6 component (e.g. a simulation, model, or package) and specifies zero or more variables.

Each variable is described by a set of attributes, some of which are required, some optional.

DFN files must contain only ASCII characters. DFN files may include comment lines, which must begin with "#".

Variables are organized into **blocks**. A block is a set of related input parameters. Each block in a DFN file is conventionally preceded by a comment line according to the format:

```
# --------------------- component subcomponent block ---------------------
```

A variable is specified by a set of attributes, with one attribute per line, following the general format:

```
name value
```

Some attribute values are optional. Variables are delimited by empty lines.

DFN files should be named `component-subcomponent.dfn`, where the component and subcomponent are abbreviated acronyms or names. For instance, a groundwater flow (GWF) model's initial conditions (IC) package DFN file is named `gwf-ic.dfn`.

DFN files are parsed to generate both source code and documentation.

**Note**: The DFN format is undergoing a migration to TOML. This involves transforming the flat variable representation into a tree with composite variables as nodes and scalars as leaves.

## Variable attributes

A MODFLOW 6 input variable is described by a set of attributes. Some attributes are optional and do not need to be specified (because they have a default value) and others are required. The following attributes are supported:

| Attribute     | Description                                                            | Required | Default | Notes                                                                                                                                                         |
|---------------|------------------------------------------------------------------------|----------|---------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| block         | Block name                                                             | Yes      |         | Each parameter must belong to a block.                                                                                                                        |
| name          | Parameter name                                                         | Yes      |         | Parameter names should be lower case and may contain (unescaped) underscores.                                                                                 |
| type          | Parameter type                                                         | Yes      |         | Valid values are: `keyword`, `string`, `integer`, `double precision`, `recarray`, `record`, `recordrepeating`, and `keystring`.                               |
| valid         | Valid parameter values                                                 | No       | None    |                                                                                                                                                               |
| shape         | Array shape or string length                                                            | No       | None    | Only required/relevant for strings or array parameters. A single value without parentheses is interpreted as string length. One or more comma-separated values enclosed by parentheses indicate array shape. Each value may a) reference another scalar variable by name (referenced variables may or may not be in the same component &mdash; e.g., it is common to reference grid dimensions defined in the DIS[V/U] package, or time dimensions defined in the TDIS package), b) reference another array variable by name, in which case the referencing variable is interpreted as the same shape as the referenced, c) reference a named constant defined in the `Constants.f90` module, or d) decline to specify a shape and specify only dimensions, with a value of `any1d`, `any2d`, or `any3d`. Values may be preceded by `<` to place an upper bound on the size of a certain dimension rather than an exact size requirement.
| tagged        | Whether a keyword is required before the parameter value.              | No       | True    | Set to false for keyword parameters (which do not take a value).                                                                                              |
| in_record     | Whether the parameter is part of a record.                             | No       | False   | If true, the parameter must follow a record parameter keyword rather than being listed on its own line.                                                       |
| layered       | Whether the parameter is layered.                                      | No       | False   | If true, then the LAYERED keyword will be written to the input instructions.                                                                                  |
| reader        | The MODFLOW 6 routine used to read the parameter value.                | No       | (inferred) | This attribute is deprecated and is now inferred from `type` and `shape` but remains supported for compatibility. Valid values are: `urword`, `readarray`. |
| optional      | Whether the parameter is optional.                                     | No       | False   |                                                                                                                                                               |
| longname      | A brief but descriptive label for the parameter.                       | Yes      |         | May contain spaces.                                                                                                                                           |
| description   | A full description of the parameter.                                   | Yes      |         | Should describe the parameter in detail. Underscores must be escaped since this value is parsed and substituted into LaTeX files for the MF6IO documentation. |
| preserve_case | Whether the case of the parameter value (if text) should be preserved. | No       | False   | This should be set to true for filename parameters.                                                                                                           |
| default_value | The parameter's default value.                                         | No       | None    | Should be a valid Python literal.                                                                                                                             |
| numeric_index | Indicates that this is an index variable.                              | No       | False   | Indicates that FloPy should treat the parameter as zero-based.                                                                                                |
| deprecated    | Indicates that the parameter has been deprecated.                      | No       | None    | Should be a semantic version number, the version in which the parameter was deprecated. If this attribute is provided without a value, it is ignored.            |
| removed       | Indicates that the parameter has been removed.                         | No       | None    | Should be a semantic version number, the version in which the parameter was removed. If this attribute is provided without a value, it is ignored.               |
| developmode    | Indicates that the parameter should be omitted from releases.               | No       | False   | Should be set to true to indicate that a parameter is not ready for inclusion in releases. The parameter will be omitted from the generated IO guide.   |

### Reader Attribute

The reader attribute indicates which read routine is used by MODFLOW 6. It is only used to generate documentation. **This attribute is now automatically inferred** from attributes `type` and `shape` and does not need to be specified in DFN files. If present in a DFN file, it will be ignored in favor of the inferred value. Array types (`integer` or `double precision` with a non-empty shape) use reader `readarray`. All other types use `urword`.

The `readarray` reader is used for array variables and results in specialized documentation.  It allows for a LAYERED keyword to be specified.  For example, the delr array of the DIS package might be specified as:

```
block griddata
name delr
type double precision
shape (ncol)
longname spacing along a row
description is the is the column spacing in the row direction.
```

This results in the following block description:

```
BEGIN GRIDDATA
  DELR
    <delr(ncol)> -- READARRAY
END GRIDDATA
```

The icelltype variable also uses readarray and is specified as:

```
block GRIDDATA
name icelltype
type integer
shape (nodes)
valid
reader readarray
layered true
optional
longname confined or convertible indicator
description flag for each cell that specifies how saturated thickness is treated.  0 means saturated thickness is held constant;  $>$0 means saturated thickness varies with computed head when head is below the cell top; $<$0 means saturated thickness varies with computed head unless the THICKSTRT option is in effect.  When THICKSTRT is in effect, a negative value of icelltype indicates that saturated thickness will be computed as STRT-BOT and held constant.
```

This results in the following block description:

```
BEGIN GRIDDATA
  ICELLTYPE [LAYERED]
    icelltype(nodes) -- READARRAY
END GRIDDATA
```

## Variable types

The following input parameter types are supported.

### Keyword

As shown in the simple example above, a keyword type is simply one that is specified using a text string.  The keyword name is automatically converted into upper case.  Upper case is used in the input instructions to indicate a word that is recognized by the MODFLOW 6 program.

A simple example of a keyword is:

```
# --------------------- gwf dis options ---------------------

block options
name nogrb
type keyword
reader urword
optional true
longname do not write binary grid file
description keyword to deactivate writing of the binary grid file.
```

which renders as,

```
BEGIN OPTIONS
  [NOGRB]
END OPTIONS
```

Note that the optional tag being set to true results in square brackets around NOGRB.

### String

As shown in the simple example above, a string type is simply a keyword followed by a text string.  This type allows the user to provide text information to MODFLOW 6.

A simple example of a keyword is:

```
# --------------------- gwf dis options ---------------------

block options
name length_units
type string
reader urword
optional true
longname model length units
description is the length units used for this model.  Values can be ``FEET'', ``METERS'', or ``CENTIMETERS''.  If not specified, the default is ``UNKNOWN''.
```

which renders as,

```
BEGIN OPTIONS
  [LENGTH_UNITS length_units]
END OPTIONS
```

Again the optional tag is used here to indicate that the keyword and string value are optional.

Though it wouldn't make sense in this context, we could also have set tagged to False.  By doing so, we would have ended up with:

```
BEGIN OPTIONS
  [length_units]
END OPTIONS
```

In this case, the user would provide length_units directly without a preceding keyword.  Eliminating tags is used later to construct records consisting of multiple data entries.

### Numeric

As you might expect, the integer and double precision types allow for specification of integer and double precision input by the user.  A simple example from the Solution Sparse Matrix Solver Nonlinear block is:

```
# --------------------- sln sms nonlinear ---------------------

block nonlinear
name outer_hclose
type double precision
reader urword
optional false
longname head change criterion
description real value defining the head change criterion for convergence of the outer (nonlinear) iterations, in units of length. When the maximum absolute value of the head change at all nodes during an iteration is less than or equal to \texttt{outer\_hclose}, iteration stops. Commonly, \texttt{outer\_hclose} equals 0.01.

block nonlinear
name outer_maximum
type integer
reader urword
optional false
longname outer maximum iterations
description integer value defining the maximum number of outer (nonlinear) iterations -- that is, calls to the solution routine. For a linear problem \texttt{outer\_maximum} should be 1.
```

This of course renders as:

```
BEGIN NONLINEAR
  OUTER_HCLOSE outer_hclose
  OUTER_MAXIMUM outer_maximum
END NONLINEAR
```

Note that in this case, both of the variables are required (optional false) and so they are not enclosed in brackets.

### Recarray

The recarray type is patterned after the recarray type that is available in the numpy package for Python.

An example of a recarray record is shown below for the drain package.  First you'll note that the recarray has a shape.  This shape is of maxbound, which is the maximum number of records that the user can enter.  Also note that following the ``recarray'' identifier is cellid, elev, cond, aux, and boundname.  These are all additional variables that are described after the recarray.  Because these are listed next to recarray, the protocol is that they will all be listed on one line.  You'll also note that the cellid, elev, cond, aux, and boundname variables have the in_record attribute set to true.  This is required so that the variables are not written again after the recarray; they are only written inside the recarray.  These variables all have the tagged attribute set to false so that they are not preceded by a keyword.

```
block period
name periodrecarray
type recarray cellid elev cond aux boundname
shape (maxbound)
reader urword
longname
description

block period
name cellid
type integer
shape (ncelldim)
tagged false
in_record true
reader urword
longname cell identifier
description REPLACE cellid {}

block period
name elev
type double precision
shape
tagged false
in_record true
reader urword
longname drain elevation
description is the elevation of the drain. If the Options block includes a TIMESERIESFILE entry (see the ``Time-Variable Input'' section), values can be obtained from a time series by entering the time-series name in place of a numeric value.

block period
name cond
type double precision
shape
tagged false
in_record true
reader urword
longname drain conductance
description is the hydraulic conductance of the interface between the aquifer and the drain. If the Options block includes a TIMESERIESFILE entry (see the ``Time-Variable Input'' section), values can be obtained from a time series by entering the time-series name in place of a numeric value.

block period
name aux
type double precision
in_record true
tagged false
shape (naux)
reader urword
optional true
longname auxiliary variables
description REPLACE aux {'{#1}': 'drain'}

block period
name boundname
type string
shape
tagged false
in_record true
reader urword
optional true
longname drain name
description REPLACE boundname {'{#1}': 'drain'}
```

When these lines are rendered, we get the following:

```
BEGIN PERIOD
  cellid(ncelldim) elev cond [aux(naux)] [boundname]
  cellid(ncelldim) elev cond [aux(naux)] [boundname]
  ...
END PERIOD
```

### Record

A record type is similar to the recarray type, except that it doesn't have a shape.  A record can be just a list of values on one line.  A simple example of this is in the SMS Linear block, where we have the rcloserecord consisting of the inner_rclose and rclose_option variables.

```
block linear
name rcloserecord
type record inner_rclose rclose_option
reader urword
optional false
longname rclose record
description

block linear
name inner_rclose
type double precision
in_record true
reader urword
optional false
longname flow residual tolerance
description real value that defines the flow residual tolerance for convergence of the SMS linear solver and specific flow residual criteria used. This value represents the maximum allowable residual at any single node.  Value is in units of length cubed per time, and must be consistent with \mf length and time units. Usually a value of $1.0 \times 10^{-1}$ is sufficient for the flow-residual criteria when meters and seconds are the defined \mf length and time.

block linear
name rclose_option
type string
tagged false
in_record true
reader urword
optional true
longname flow residual tolerance
description an optional keyword that defines the specific flow residual criterion used.  \texttt{L2NORM\_RCLOSE}--an optional keyword that is used to specify that \texttt{inner\_rclose} represents a L-2 Norm closure criteria instead of a infinity-Norm (absolute convergence criteria). When \texttt{L2NORM\_RCLOSE} is specified, a reasonable initial \texttt{inner\_rclose} value is $\left( 1.0 \times 10^{-1} \times \text{active nodes} \right)$ when meters and seconds are the defined \mf length and time.  \texttt{RELATIVE\_RCLOSE}--an optional keyword that is used to specify that \texttt{inner\_rclose} represents a relative L-2 Norm reduction closure criteria instead of a infinity-Norm (absolute convergence criteria). When \texttt{RELATIVE\_RCLOSE} is specified, a reasonable initial \texttt{inner\_rclose} value is $1.0 \times 10^{-4}$ and convergence is achieved for a given inner (linear) iteration when $\Delta h \le$ \texttt{inner\_hclose} and the current L-2 Norm is $\le$ the product of the \texttt{RELATIVE\_RCLOSE} and the initial L-2 Norm for the current inner (linear) iteration. If \texttt{rclose\_option} is not specified, an absolute residual (infinity-norm) criterion is used.
```

This renders as the following.  The inner_rclose variable is tagged by default and so the INNER_RCLOSE keyword precedes inner_rclose.  The rclose_option variable has tagged set to false, and so there is no keyword preceding it.  It also has optional set to true, and so there are square brackets around it.

```
BEGIN LINEAR
  INNER_RCLOSE inner_rclose [rclose_option]
END LINEAR
```

#### Repeating records

The recordrepeating type is the same as the record type.  We've given it another name to indicate that the user can specify more than on line.  The only place this is used at present is in the Output Control.


### Keystring

A keystring is a union type.  A keystring identifies that a variable can specified using multiple different types of input.  It is most often used by the advanced packages to adjust individual comments of the package.  The following is an example of the mawsetting variable:

```
block period
name mawsetting
type keystring status flowing_wellrecord rate well_head head_limit shutoffrecord rate_scalingrecord auxiliaryrecord
shape
tagged false
in_record true
reader urword
longname
description line of information that is parsed into a keyword and values.  Keyword values that can be used to start the \texttt{mawsetting} string include: \texttt{FLOWING\_WELL}, \texttt{RATE}, \texttt{WELL\_HEAD}, \texttt{HEAD\_LIMIT}, \texttt{SHUT\_OFF}, \texttt{RATE\_SCALING}, and \texttt{AUXILIARY}.

block period
name status
type string
shape
tagged true
in_record true
reader urword
longname well status
description keyword option to define well status.  \texttt{status} can be \texttt{ACTIVE}, \texttt{INACTIVE}, or \texttt{CONSTANT}. By default, \texttt{status} is \texttt{ACTIVE}.

block period
name flowing_wellrecord
type record flowing_well fwelev fwcond
shape
tagged
in_record true
reader urword
longname
description

block period
name flowing_well
type keyword
shape
in_record true
reader urword
longname well is a flowing well
description keyword to indicate the well is a flowing well.  The \texttt{flowing\_well} option can be used to simulate flowing wells when the simulated well head exceeds the specified drainage elevation.

block period
name fwelev
type double precision
shape
tagged false
in_record true
reader urword
longname flowing well elevation
description elevation used to determine whether or not the well is flowing.

block period
name fwcond
type double precision
shape
tagged false
in_record true
reader urword
longname well flowing well elevation
description elevation used to determine whether or not the well is flowing.

block period
name rate
type double precision
shape
tagged true
in_record true
reader urword
longname well pumping rate
description is the volumetric pumping rate for the multi-aquifer well. A positive value indicates recharge and a negative value indicates discharge (pumping). \texttt{rate} only applies to active (\texttt{IBOUND}$>0$) multi-aquifer wells. If the Options block includes a \texttt{TIMESERIESFILE} entry (see the ``Time-Variable Input'' section), values can be obtained from a time series by entering the time-series name in place of a numeric value. By default, the \texttt{rate} for each multi-aquifer well is zero.

block period
name well_head
type double precision
shape
tagged true
in_record true
reader urword
longname well head
description is the head in the multi-aquifer well. \texttt{well\_head} is only applied to constant head (\texttt{STATUS} is \texttt{CONSTANT}) and inactive (\texttt{STATUS} is \texttt{INACTIVE}) multi-aquifer wells.

block period
name head_limit
type string
shape
tagged true
in_record true
reader urword
longname head limit
description is the limiting water level (head) in the well, which is the minimum of the well \texttt{rate} or the well inflow rate from the aquifer. \texttt{head\_limit} is only applied to discharging wells (\texttt{rate}$<0$). \texttt{head\_limit} can be deactivated by specifying the text string `\texttt{off}'. By default, \texttt{head\_limit} is `\texttt{off}'.

block period
name shutoffrecord
type record shut_off minrate maxrate
shape
tagged
in_record true
reader urword
longname
description

block period
name shut_off
type keyword
shape
in_record true
reader urword
longname shut off well
description keyword for activating well shut off capability.  Subsequent values define the minimum and maximum pumping rate that a well must exceed to shutoff or reactivate a well, respectively, during a stress period. \texttt{shut\_off} is only applied to discharging wells (\texttt{rate}$<0$) and if \texttt{head\_limit} is specified (not set to `\texttt{off}').  If \texttt{head\_limit} is specified, \texttt{shut\_off} can be deactivated by specifying a minimum value equal to zero. By default, \texttt{shut\_off} is not used.

block period
name minrate
type double precision
shape
tagged false
in_record true
reader urword
longname minimum shutoff rate
description is the minimum rate that a well must exceed to shutoff a well during a stress period.  \texttt{minrate} must be less than \texttt{maxrate}.

block period
name maxrate
type double precision
shape
tagged false
in_record true
reader urword
longname maximum shutoff rate
description is the maximum rate that a well must exceed to reactivate a well during a stress period.  \texttt{maxrate} must be less than \texttt{minrate}.

block period
name rate_scalingrecord
type record rate_scaling pump_elevation scaling_length
shape
tagged
in_record true
reader urword
longname
description

block period
name rate_scaling
type keyword
shape
in_record true
reader urword
longname rate scaling
description activate rate scaling.  If \texttt{rate\_scaling} is specified, both \texttt{pump\_elevation} and \texttt{reduction\_length} must be specified. \texttt{rate\_scaling} cannot be used with \texttt{head\_limit}.

block period
name pump_elevation
type double precision
shape
tagged false
in_record true
reader urword
longname pump elevation
description is the elevation of the multi-aquifer well pump (\texttt{pump\_elevation}).  \texttt{pump\_elevation} cannot be less than the bottom elevation (\texttt{bottom}) of the multi-aquifer well. By default, \texttt{pump\_elevation} is set equal to the bottom of the largest \texttt{GWF} node number connected to a MAW well.

block period
name scaling_length
type double precision
shape
tagged false
in_record true
reader urword
longname
description height above the pump elevation (\texttt{scaling\_length}) below which the pumping rate is reduced.  The default value for \texttt{scaling\_length} is the well radius.

block period
name auxiliaryrecord
type record auxiliary auxname auxval
shape
tagged
in_record true
reader urword
longname
description

block period
name auxiliary
type keyword
shape
in_record true
reader urword
longname
description keyword for specifying auxiliary variable.

block period
name auxname
type string
shape
tagged false
in_record true
reader urword
longname
description name for the auxiliary variable to be assigned \texttt{auxval}.  \texttt{auxname} must match one of the auxiliary variable names defined in the \texttt{OPTIONS} block. If \texttt{auxname} does not match one of the auxiliary variable names defined in the \texttt{OPTIONS} block the data are ignored.

block period
name auxval
type double precision
shape
tagged false
in_record true
reader urword
longname auxiliary variable value
description value for the auxiliary variable.
```

The mawsetting variable renders very simply in the block description as

```
BEGIN PERIOD iper
  WELL wellno mawsetting
  WELL wellno mawsetting
  ...
END PERIOD
```

But, the description for the mawsetting contains the additional information specified for the keystring.  For example, the following is written as the latex description for mawsetting:

```
\item \texttt{mawsetting}---line of information that is parsed into a keyword and values.  Keyword values that can be used to start the \texttt{mawsetting} string include: \texttt{FLOWING\_WELL}, \texttt{RATE}, \texttt{WELL\_HEAD}, \texttt{HEAD\_LIMIT}, \texttt{SHUT\_OFF}, \texttt{RATE\_SCALING}, and \texttt{AUXILIARY}.

\begin{verbatim}
STATUS status
FLOWING_WELL fwelev fwcond
RATE rate
WELL_HEAD well_head
HEAD_LIMIT head_limit
SHUT_OFF minrate maxrate
RATE_SCALING pump_elevation scaling_length
AUXILIARY auxname auxval
\end{verbatim}
```



## Text substitutions

Many of the MODFLOW 6 input variables share common description information.  This information can be defined in one place and then referenced as many times as needed throughout the definition files.  The definition file [./dfn/common.dfn](./dfn/common.dfn) is where common information is defined.  For example, there is a text string defined as:

```
name auxnames
description defines an array of one or more auxiliary variable names.  There is no limit on the number of auxiliary variables that can be provided on this line; however, lists of information provided in subsequent blocks must have a column of data for each auxiliary variable name defined here.   The number of auxiliary variables detected on this line determines the value for naux.  Comments cannot be provided anywhere on this line as they will be interpreted as auxiliary variable names.  Auxiliary variables may not be used by the {#1} model, but they will be available for use by other parts of the program.  The ``AUX'' and ``AUXILIARY'' keywords can be used as a substitute for ``AUXNAMES''. The program will terminate with an error if auxiliary variables are specified on more than one line in the options block.
```

This can be used in the description for auxnames throughout in the following manner:

```
block options
name auxnames
type string
shape (naux)
reader urword
optional true
longname keyword to save GWFGWF flows
description REPLACE auxnames {'{#1}': 'Groundwater Flow'}
```

In the description attribute, the capital REPLACE instructs the processor to replace auxnames with the text string defined by auxnames in common.dfn.  Also included here is a Python dictionary, which instructs the processor to replace the text string '{#1}' with 'Groundwater Flow'.

## Writing definition files

This section demonstrates a definition file for a MODFLOW 6 component.

Below is shown DFN file contents specifying an options block with optional keywords:

```
# --------------------- gwf dis options ---------------------

block options
name length_units
type string
reader urword
optional true
longname model length units
description is the length units used for this model.  Values can be ``FEET'', ``METERS'', or ``CENTIMETERS''.  If not specified, the default is ``UNKNOWN''.

block options
name nogrb
type keyword
reader urword
optional true
longname do not write binary grid file
description keyword to deactivate writing of the binary grid file.
```

This corresponds to the following in a MODFLOW 6 input file:

```
BEGIN OPTIONS
  [LENGTH_UNITS length_units]
  [NOGRB]
END OPTIONS
```

Square brackets indicate that a parameter is optional.

The MODFLOW 6 repository contains scripts to generate the above snippet the DFN file and substitute it into the MODFLOW 6 input/output documentation. In addition, a LaTeX file will be created, the name of which will be the prefix of the definition file plus `-desc.tex`.  This LaTeX file describes each parameter in the definition file:

```
\item \texttt{length\_units}---is the length units used for this model.  Values can be ``FEET'', ``METERS'', or ``CENTIMETERS''.  If not specified, the default is ``UNKNOWN''.

\item \texttt{NOGRB}---keyword to deactivate writing of the binary grid file.
```

The parameters will also be listed in a comprehensive Markdown table, which is also inserted into the MF6IO documentation and hosted on the online documentation site:


| component | package | block | variable name | type | description |
| :---: | :---: | :---: | :---: | :---: | --- |
| GWF | DIS | OPTIONS | LENGTH_UNITS | STRING | is the length units used for this model.  Values can be ``FEET'', ``METERS'', or ``CENTIMETERS''.  If not specified, the default is ``UNKNOWN''. |
| GWF | DIS | OPTIONS | NOGRB | KEYWORD | keyword to deactivate writing of the binary grid file. |

Further consideration is needed when writing definition files for certain MODFLOW 6 components.

### Subpackages

A subpackage is a package referenced by another package (vs being referenced in the name file).  The tas, ts, and obs packages are examples of subpackages.  There are a few additional steps required when creating a subpackage definition file.  First, verify that the parent package's dfn file has a file record for the subpackage to the option block.   For example, for the time series package the file record definition starts with:

```
block options
name ts_filerecord
type record ts6 filein ts6_filename
```

Verify that the same naming convention is followed as the example above, specifically:

```
name <subpackage-abbr>_filerecord
record <subpackage-abbr>6 filein <subpackage-abbr>6_filename
```

Next, create the child package definition file (see sections below for more information). 

When your child package is ready for release follow the same procedure described in  "Creating a Definition File that will Work with Flopy" above along with a few additional steps required for subpackages.  At the top of the child dfn file add two lines describing how the parent and child packages are related. The first line determines how the subpackage is linked to the package:

```
# flopy subpackage \<parent record\> \<abbreviation\> \<child data\> \<data name\>
```

* Parent record is the MF6 record name of the filerecord in parent package that references the child packages file name
* Abbreviation is the short abbreviation of the new subclass
* Child data is the name of the child class data that can be passed in as parameter to the parent class. Passing in this parameter to the parent class automatically creates the child class with the data provided.
* Data name is the parent class parameter name that automatically creates the child class with the data provided.

The example below is the first line from the ts subpackage dfn:

```
# flopy subpackage ts_filerecord ts timeseries timeseries
```

The second line determines the variable name of the subpackage's parent and the type of parent (the parent package's object oriented parent):

```
# flopy parent_name_type \<parent package variable name\> \<parent package type\>
```

An example below is the second line in the ts subpackage dfn:

```
# flopy parent_name_type parent_package MFPackage
```

There are three possible types (or combination of them) that can be used for "parent package type", MFPackage, MFModel, and MFSimulation. If a package supports multiple types of parents (for example, it can be either in the model namefile or in a package, like the obs package), include all the types supported, separating each type with a / (MFPackage/MFModel).

### Models

To create a new type of model choose a unique three letter model abbreviation ("gwf", "gwt", ...). Create a name file dfn with the naming convention \<model abbr\>-nam.dfn. The name file must have only an options and packages block (see gwf-nam.dfn as an example). Create a new dfn file for each of the packages in your new model, following the naming convention described above. 

When your model is ready for release copy the dfn file to the flopy distribution in the flopy/mf6/data/dfn folder, run createpackages.py, and check in your new dfn files, the package classes, and updated init.py that createpackages.py created.

### Solvers

Create a solver definition file as you would any package definition file.  When you are done add a commented line at the top of the definition file to let FloPy know that this package is a solver (solution package type). The line should look like this:

```
# solution_package <solver abbreviation> <list of model abbreviations the solver supports>
```

For example, the following would tell FloPy the IMS package supports the gwf6 and gwt6 model types:

```
# flopy solution_package ims gwf6 gwt6
```

If a "*" is used instead of the list of model abbreviations, the solver is assumed to support all model types.  For example, the following would tell FloPy that the IMS package supports all models:

```
# flopy solution_package ims *
```

## Developer tooling

As mentioned, the MODFLOW 6 repository contains scripts to generate source code and documentation from DFN files.

The Python script `mf6ivar.py` will process all of the definition files and create a Markdown file, LaTeX files of the variable descriptions, and text files containing the blocks.

The Python script `deprecations.py` will search definition files for `deprecated` or `removed` options and create a Markdown file containing a table of deprecations and removals. The Python script `mk_deprecations.py` will convert that table to LaTeX for the PDF documentation.

The Python script `dfn2f90.py` will generate a Fortran interface layer for the IDM-integrated components.
! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfLakInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_lak_param_definitions
  public gwf_lak_aggregate_definitions
  public gwf_lak_block_definitions
  public GwfLakParamFoundType
  public gwf_lak_multi_package
  public gwf_lak_subpackages

  type GwfLakParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: print_stage = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: stage_filerecord = .false.
    logical :: options_stage = .false.
    logical :: stagefile = .false.
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: pkgconvfilerec = .false.
    logical :: pkgconv = .false.
    logical :: pkgconvfname = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: options_filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: iforcefb = .false.
    logical :: igwhcopt = .false.
    logical :: outdmax = .false.
    logical :: surfdep = .false.
    logical :: iimplicit = .false.
    logical :: maxlakit = .false.
    logical :: dmaxchg = .false.
    logical :: convtime = .false.
    logical :: convlength = .false.
    logical :: iconvchk = .false.
    logical :: nlakes = .false.
    logical :: noutlets = .false.
    logical :: ntables = .false.
    logical :: packagedata_ifno = .false.
    logical :: strt = .false.
    logical :: nlakeconn = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: conndata_ifno = .false.
    logical :: iconn = .false.
    logical :: cellid = .false.
    logical :: claktype = .false.
    logical :: bedleak = .false.
    logical :: belev = .false.
    logical :: telev = .false.
    logical :: connlen = .false.
    logical :: connwidth = .false.
    logical :: tables_ifno = .false.
    logical :: tab6 = .false.
    logical :: tables_filein = .false.
    logical :: tab6_filename = .false.
    logical :: outlets_outletno = .false.
    logical :: lakein = .false.
    logical :: lakeout = .false.
    logical :: couttype = .false.
    logical :: outlets_invert = .false.
    logical :: outlets_width = .false.
    logical :: outlets_rough = .false.
    logical :: outlets_slope = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: stage_in = .false.
    logical :: rainfall_in = .false.
    logical :: evaporation_in = .false.
    logical :: runoff_in = .false.
    logical :: inflow_in = .false.
    logical :: withdrawal_in = .false.
    logical :: rate_in = .false.
    logical :: invert_in = .false.
    logical :: rough_in = .false.
    logical :: width_in = .false.
    logical :: slope_in = .false.
    logical :: auxrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwfLakParamFoundType

  logical :: gwf_lak_multi_package = .true.

  character(len=16), parameter :: &
    gwf_lak_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwflak_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_print_stage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_STAGE', & ! tag name
    'PRINT_STAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated stages to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save lake flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_stage_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE_FILERECORD', & ! tag name
    'STAGE_FILERECORD', & ! fortran variable
    'RECORD STAGE FILEOUT STAGEFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_options_stage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE', & ! tag name
    'OPTIONS_STAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'stage keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_stagefile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'STAGEFILE', & ! tag name
    'STAGEFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET_FILERECORD', & ! tag name
    'BUDFILEREC', & ! fortran variable
    'RECORD BUDGET FILEOUT BUDGETFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET', & ! tag name
    'BUDGET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETFILE', & ! tag name
    'BUDGETFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV_FILERECORD', & ! tag name
    'BUDCSVFILEREC', & ! fortran variable
    'RECORD BUDGETCSV FILEOUT BUDGETCSVFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV', & ! tag name
    'BUDGETCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSVFILE', & ! tag name
    'BUDGETCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_pkgconvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILERECORD', & ! tag name
    'PKGCONVFILEREC', & ! fortran variable
    'RECORD PACKAGE_CONVERGENCE FILEOUT '// &
    'PACKAGE_CONVERGENCE_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_pkgconv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE', & ! tag name
    'PKGCONV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'package_convergence keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_pkgconvfname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILENAME', & ! tag name
    'PKGCONVFNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_options_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'OPTIONS_FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'MOVER', & ! tag name
    'MOVER', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iforcefb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_FORCE_FALLBACK', & ! tag name
    'IFORCEFB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'force the substitution fallback solution', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_igwhcopt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_GROUNDWATER_HEAD_CONDUCTANCE', & ! tag name
    'IGWHCOPT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use the groundwater head to calculate horizontal conductance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outdmax = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_MAXIMUM_OUTLET_DEPTH', & ! tag name
    'OUTDMAX', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'maximum outlet depth', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_surfdep = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'SURFDEP', & ! tag name
    'SURFDEP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'surface depression depth', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iimplicit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'IMPLICIT', & ! tag name
    'IIMPLICIT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'solve lake stage in the groundwater flow matrix', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_maxlakit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'MAXIMUM_ITERATIONS', & ! tag name
    'MAXLAKIT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'LAK Newton-Raphson iterations', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_dmaxchg = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'MAXIMUM_STAGE_CHANGE', & ! tag name
    'DMAXCHG', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'stage closure tolerance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_convtime = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_CONVERSION', & ! tag name
    'CONVTIME', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'time conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_convlength = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_CONVERSION', & ! tag name
    'CONVLENGTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'length conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iconvchk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_FINAL_CHECK', & ! tag name
    'ICONVCHK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'disable final convergence check for lake stages', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_nlakes = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'DIMENSIONS', & ! block
    'NLAKES', & ! tag name
    'NLAKES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of lakes', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_noutlets = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'DIMENSIONS', & ! block
    'NOUTLETS', & ! tag name
    'NOUTLETS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of outlets', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ntables = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'DIMENSIONS', & ! block
    'NTABLES', & ! tag name
    'NTABLES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of tables', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_strt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting lake stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_nlakeconn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NLAKECONN', & ! tag name
    'NLAKECONN', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of lake connections', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_aux = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'AUX', & ! tag name
    'AUX', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_conndata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IFNO', & ! tag name
    'CONNDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_iconn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'ICONN', & ! tag name
    'ICONN', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'connection number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_claktype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CLAKTYPE', & ! tag name
    'CLAKTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake connection type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_bedleak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'BEDLEAK', & ! tag name
    'BEDLEAK', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'bed leakance', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_belev = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'BELEV', & ! tag name
    'BELEV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_telev = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'TELEV', & ! tag name
    'TELEV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'top elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_connlen = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CONNLEN', & ! tag name
    'CONNLEN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection length', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_connwidth = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CONNWIDTH', & ! tag name
    'CONNWIDTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_tables_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'TABLES', & ! block
    'IFNO', & ! tag name
    'TABLES_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_tab6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'TABLES', & ! block
    'TAB6', & ! tag name
    'TAB6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_tables_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'TABLES', & ! block
    'FILEIN', & ! tag name
    'TABLES_FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_tab6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'TABLES', & ! block
    'TAB6_FILENAME', & ! tag name
    'TAB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'table file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets_outletno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'OUTLETNO', & ! tag name
    'OUTLETS_OUTLETNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'outlet number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_lakein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'LAKEIN', & ! tag name
    'LAKEIN', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for upstream lake', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_lakeout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'LAKEOUT', & ! tag name
    'LAKEOUT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for downstream lake', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_couttype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'COUTTYPE', & ! tag name
    'COUTTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'outlet type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets_invert = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'INVERT', & ! tag name
    'OUTLETS_INVERT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'invert elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets_width = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'WIDTH', & ! tag name
    'OUTLETS_WIDTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'outlet width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets_rough = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'ROUGH', & ! tag name
    'OUTLETS_ROUGH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets_slope = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'SLOPE', & ! tag name
    'OUTLETS_SLOPE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'bed slope', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'NUMBER', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake or outlet number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_status = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_stage_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'STAGE', & ! tag name
    'STAGE_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'lake stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_rainfall_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'RAINFALL', & ! tag name
    'RAINFALL_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'rainfall rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_evaporation_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'EVAPORATION', & ! tag name
    'EVAPORATION_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'evaporation rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_runoff_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'RUNOFF', & ! tag name
    'RUNOFF_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'runoff rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_inflow_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'INFLOW', & ! tag name
    'INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'inflow rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_withdrawal_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'WITHDRAWAL', & ! tag name
    'WITHDRAWAL_IN', & ! fortran variable
    'STRING', & ! type
    'NLAKES', & ! shape
    'maximum withdrawal rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_rate_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE_IN', & ! fortran variable
    'STRING', & ! type
    'NOUTLETS', & ! shape
    'extraction rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_invert_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'INVERT', & ! tag name
    'INVERT_IN', & ! fortran variable
    'STRING', & ! type
    'NOUTLETS', & ! shape
    'invert elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_rough_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'ROUGH', & ! tag name
    'ROUGH_IN', & ! fortran variable
    'STRING', & ! type
    'NOUTLETS', & ! shape
    'roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_width_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'WIDTH', & ! tag name
    'WIDTH_IN', & ! fortran variable
    'STRING', & ! type
    'NOUTLETS', & ! shape
    'outlet width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_slope_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'SLOPE', & ! tag name
    'SLOPE_IN', & ! fortran variable
    'STRING', & ! type
    'NOUTLETS', & ! shape
    'bed slope', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_auxrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARYRECORD', & ! tag name
    'AUXRECORD', & ! fortran variable
    'RECORD AUXILIARY AUXNAME AUXVAL', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARY', & ! tag name
    'PERIOD_AUXILIARY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_auxname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'AUXNAME', & ! tag name
    'AUXNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_auxval = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'AUXVAL', & ! tag name
    'AUXVAL', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'auxiliary variable value', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_lak_param_definitions(*) = &
    [ &
    gwflak_auxiliary, &
    gwflak_boundnames, &
    gwflak_iprpak, &
    gwflak_print_stage, &
    gwflak_iprflow, &
    gwflak_ipakcb, &
    gwflak_stage_filerecord, &
    gwflak_options_stage, &
    gwflak_stagefile, &
    gwflak_budfilerec, &
    gwflak_budget, &
    gwflak_fileout, &
    gwflak_budgetfile, &
    gwflak_budcsvfilerec, &
    gwflak_budgetcsv, &
    gwflak_budgetcsvfile, &
    gwflak_pkgconvfilerec, &
    gwflak_pkgconv, &
    gwflak_pkgconvfname, &
    gwflak_ts_filerecord, &
    gwflak_ts6, &
    gwflak_options_filein, &
    gwflak_ts6_filename, &
    gwflak_obs_filerecord, &
    gwflak_obs6, &
    gwflak_obs6_filename, &
    gwflak_mover, &
    gwflak_iforcefb, &
    gwflak_igwhcopt, &
    gwflak_outdmax, &
    gwflak_surfdep, &
    gwflak_iimplicit, &
    gwflak_maxlakit, &
    gwflak_dmaxchg, &
    gwflak_convtime, &
    gwflak_convlength, &
    gwflak_iconvchk, &
    gwflak_nlakes, &
    gwflak_noutlets, &
    gwflak_ntables, &
    gwflak_packagedata_ifno, &
    gwflak_strt, &
    gwflak_nlakeconn, &
    gwflak_aux, &
    gwflak_boundname, &
    gwflak_conndata_ifno, &
    gwflak_iconn, &
    gwflak_cellid, &
    gwflak_claktype, &
    gwflak_bedleak, &
    gwflak_belev, &
    gwflak_telev, &
    gwflak_connlen, &
    gwflak_connwidth, &
    gwflak_tables_ifno, &
    gwflak_tab6, &
    gwflak_tables_filein, &
    gwflak_tab6_filename, &
    gwflak_outlets_outletno, &
    gwflak_lakein, &
    gwflak_lakeout, &
    gwflak_couttype, &
    gwflak_outlets_invert, &
    gwflak_outlets_width, &
    gwflak_outlets_rough, &
    gwflak_outlets_slope, &
    gwflak_ifno, &
    gwflak_status, &
    gwflak_stage_in, &
    gwflak_rainfall_in, &
    gwflak_evaporation_in, &
    gwflak_runoff_in, &
    gwflak_inflow_in, &
    gwflak_withdrawal_in, &
    gwflak_rate_in, &
    gwflak_invert_in, &
    gwflak_rough_in, &
    gwflak_width_in, &
    gwflak_slope_in, &
    gwflak_auxrecord, &
    gwflak_period_auxiliary, &
    gwflak_auxname, &
    gwflak_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwflak_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO STRT NLAKECONN AUX BOUNDNAME', & ! type
    'NLAKES', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_connectiondata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CONNECTIONDATA', & ! tag name
    'CONNECTIONDATA', & ! fortran variable
    'RECARRAY IFNO ICONN CELLID CLAKTYPE BEDLEAK BELEV TELEV '// &
    'CONNLEN CONNWIDTH', & ! type
    'NLAKECONN', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_tables = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'TABLES', & ! block
    'TABLES', & ! tag name
    'TABLES', & ! fortran variable
    'RECARRAY IFNO TAB6 FILEIN TAB6_FILENAME', & ! type
    'NTABLES', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_outlets = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'OUTLETS', & ! block
    'OUTLETS', & ! tag name
    'OUTLETS', & ! fortran variable
    'RECARRAY OUTLETNO LAKEIN LAKEOUT COUTTYPE INVERT WIDTH ROUGH '// &
    'SLOPE', & ! type
    'NOUTLETS', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY NUMBER LAKSETTING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwflak_laksetting = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'LAK', & ! subcomponent
    'PERIOD', & ! block
    'LAKSETTING', & ! tag name
    'LAKSETTING', & ! fortran variable
    'KEYSTRING STATUS STAGE RAINFALL EVAPORATION RUNOFF INFLOW '// &
    'WITHDRAWAL RATE INVERT WIDTH SLOPE ROUGH AUXILIARYRECORD', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_lak_aggregate_definitions(*) = &
    [ &
    gwflak_packagedata, &
    gwflak_connectiondata, &
    gwflak_tables, &
    gwflak_outlets, &
    gwflak_perioddata, &
    gwflak_laksetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_lak_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CONNECTIONDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'TABLES', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'OUTLETS', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfLakInputModule

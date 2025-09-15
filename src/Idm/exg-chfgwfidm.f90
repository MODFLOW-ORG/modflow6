! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgChfgwfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_chfgwf_param_definitions
  public exg_chfgwf_aggregate_definitions
  public exg_chfgwf_block_definitions
  public ExgChfgwfParamFoundType
  public exg_chfgwf_multi_package
  public exg_chfgwf_subpackages

  type ExgChfgwfParamFoundType
    logical :: ipr_input = .false.
    logical :: ipr_flow = .false.
    logical :: ifixedcond = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: filein = .false.
    logical :: obs6_filename = .false.
    logical :: nexg = .false.
    logical :: cellidm1 = .false.
    logical :: cellidm2 = .false.
    logical :: bedleak = .false.
    logical :: cfact = .false.
  end type ExgChfgwfParamFoundType

  logical :: exg_chfgwf_multi_package = .true.

  character(len=16), parameter :: &
    exg_chfgwf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_ipr_input = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPR_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print input to list file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_ipr_flow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPR_FLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print chfgwf flows to list file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_ifixedcond = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'FIXED_CONDUCTANCE', & ! tag name
    'IFIXEDCOND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate conductance is fixed', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_obs_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_obs6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_filein = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_obs6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NEXG', & ! tag name
    'NEXG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of exchanges', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM1', & ! tag name
    'CELLIDM1', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of cell in surface water model', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM2', & ! tag name
    'CELLIDM2', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of cell in groundwater model', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_bedleak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'BEDLEAK', & ! tag name
    'BEDLEAK', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'bed leakance', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_cfact = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CFACT', & ! tag name
    'CFACT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'factor used for conductance calculation', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_chfgwf_param_definitions(*) = &
    [ &
    exgchfgwf_ipr_input, &
    exgchfgwf_ipr_flow, &
    exgchfgwf_ifixedcond, &
    exgchfgwf_obs_filerecord, &
    exgchfgwf_obs6, &
    exgchfgwf_filein, &
    exgchfgwf_obs6_filename, &
    exgchfgwf_nexg, &
    exgchfgwf_cellidm1, &
    exgchfgwf_cellidm2, &
    exgchfgwf_bedleak, &
    exgchfgwf_cfact &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgchfgwf_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'CHFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'EXCHANGEDATA', & ! tag name
    'EXCHANGEDATA', & ! fortran variable
    'RECARRAY CELLIDM1 CELLIDM2 BEDLEAK CFACT', & ! type
    'NEXG', & ! shape
    'exchange data', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_chfgwf_aggregate_definitions(*) = &
    [ &
    exgchfgwf_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_chfgwf_block_definitions(*) = &
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
    'EXCHANGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ExgChfgwfInputModule

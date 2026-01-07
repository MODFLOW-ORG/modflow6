! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgOlfgwfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_olfgwf_param_definitions
  public exg_olfgwf_aggregate_definitions
  public exg_olfgwf_block_definitions
  public ExgOlfgwfParamFoundType
  public exg_olfgwf_multi_package
  public exg_olfgwf_subpackages

  type ExgOlfgwfParamFoundType
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
  end type ExgOlfgwfParamFoundType

  logical :: exg_olfgwf_multi_package = .true.

  character(len=16), parameter :: &
    exg_olfgwf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_ipr_input = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPR_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print input to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_ipr_flow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPR_FLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print olfgwf flows to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_ifixedcond = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'FIXED_CONDUCTANCE', & ! tag name
    'IFIXEDCOND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate conductance is fixed', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_obs_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
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
    exgolfgwf_obs6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
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
    exgolfgwf_filein = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
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
    exgolfgwf_obs6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
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
    exgolfgwf_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NEXG', & ! tag name
    'NEXG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of exchanges', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM1', & ! tag name
    'CELLIDM1', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of cell in surface water model', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM2', & ! tag name
    'CELLIDM2', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of cell in groundwater model', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_bedleak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'BEDLEAK', & ! tag name
    'BEDLEAK', & ! fortran variable
    'DOUBLE', & ! type
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
    exgolfgwf_cfact = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CFACT', & ! tag name
    'CFACT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'factor used for conductance calculation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_olfgwf_param_definitions(*) = &
    [ &
    exgolfgwf_ipr_input, &
    exgolfgwf_ipr_flow, &
    exgolfgwf_ifixedcond, &
    exgolfgwf_obs_filerecord, &
    exgolfgwf_obs6, &
    exgolfgwf_filein, &
    exgolfgwf_obs6_filename, &
    exgolfgwf_nexg, &
    exgolfgwf_cellidm1, &
    exgolfgwf_cellidm2, &
    exgolfgwf_bedleak, &
    exgolfgwf_cfact &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgolfgwf_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'OLFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'EXCHANGEDATA', & ! tag name
    'EXCHANGEDATA', & ! fortran variable
    'RECARRAY CELLIDM1 CELLIDM2 BEDLEAK CFACT', & ! type
    'NEXG', & ! shape
    'exchange data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_olfgwf_aggregate_definitions(*) = &
    [ &
    exgolfgwf_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_olfgwf_block_definitions(*) = &
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

end module ExgOlfgwfInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfStoInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_sto_param_definitions
  public chf_sto_aggregate_definitions
  public chf_sto_block_definitions
  public ChfStoParamFoundType
  public chf_sto_multi_package
  public chf_sto_subpackages

  type ChfStoParamFoundType
    logical :: ipakcb = .false.
    logical :: export_ascii = .false.
    logical :: steady_state = .false.
    logical :: transient = .false.
  end type ChfStoParamFoundType

  logical :: chf_sto_multi_package = .false.

  character(len=16), parameter :: &
    chf_sto_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfsto_ipakcb = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save NPF flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfsto_export_ascii = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfsto_steady_state = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'STO', & ! subcomponent
    'PERIOD', & ! block
    'STEADY-STATE', & ! tag name
    'STEADY_STATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'steady state indicator', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfsto_transient = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'STO', & ! subcomponent
    'PERIOD', & ! block
    'TRANSIENT', & ! tag name
    'TRANSIENT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'transient indicator', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_sto_param_definitions(*) = &
    [ &
    chfsto_ipakcb, &
    chfsto_export_ascii, &
    chfsto_steady_state, &
    chfsto_transient &
    ]

  type(InputParamDefinitionType), parameter :: &
    chf_sto_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_sto_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module ChfStoInputModule

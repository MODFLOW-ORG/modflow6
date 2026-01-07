! ** Do Not Modify! MODFLOW 6 system generated file. **
module SwfStoInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public swf_sto_param_definitions
  public swf_sto_aggregate_definitions
  public swf_sto_block_definitions
  public SwfStoParamFoundType
  public swf_sto_multi_package
  public swf_sto_subpackages

  type SwfStoParamFoundType
    logical :: ipakcb = .false.
    logical :: export_ascii = .false.
    logical :: steady_state = .false.
    logical :: transient = .false.
  end type SwfStoParamFoundType

  logical :: swf_sto_multi_package = .false.

  character(len=16), parameter :: &
    swf_sto_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfsto_ipakcb = InputParamDefinitionType &
    ( &
    'SWF', & ! component
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
    swfsto_export_ascii = InputParamDefinitionType &
    ( &
    'SWF', & ! component
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
    swfsto_steady_state = InputParamDefinitionType &
    ( &
    'SWF', & ! component
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
    swfsto_transient = InputParamDefinitionType &
    ( &
    'SWF', & ! component
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
    swf_sto_param_definitions(*) = &
    [ &
    swfsto_ipakcb, &
    swfsto_export_ascii, &
    swfsto_steady_state, &
    swfsto_transient &
    ]

  type(InputParamDefinitionType), parameter :: &
    swf_sto_aggregate_definitions(*) = &
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
    swf_sto_block_definitions(*) = &
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

end module SwfStoInputModule

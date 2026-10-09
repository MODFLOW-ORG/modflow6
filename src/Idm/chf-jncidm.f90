! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfJncInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_jnc_param_definitions
  public chf_jnc_aggregate_definitions
  public chf_jnc_block_definitions
  public ChfJncParamFoundType
  public chf_jnc_multi_package
  public chf_jnc_is_advanced
  public chf_jnc_subpackages

  type ChfJncParamFoundType
    logical :: ipakcb = .false.
  end type ChfJncParamFoundType

  logical :: chf_jnc_multi_package = .false.
  logical :: chf_jnc_is_advanced = .false.

  character(len=16), parameter :: &
    chf_jnc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfjnc_ipakcb = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'JNC', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save junction flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_jnc_param_definitions(*) = &
    [ &
    chfjnc_ipakcb &
    ]

  type(InputParamDefinitionType), parameter :: &
    chf_jnc_aggregate_definitions(*) = &
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
    chf_jnc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ChfJncInputModule

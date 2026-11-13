! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgSwiswiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_swiswi_param_definitions
  public exg_swiswi_aggregate_definitions
  public exg_swiswi_block_definitions
  public ExgSwiswiParamFoundType
  public exg_swiswi_multi_package
  public exg_swiswi_is_advanced
  public exg_swiswi_subpackages

  type ExgSwiswiParamFoundType
    logical :: inocrossstorage = .false.
    logical :: inocrossflow = .false.
  end type ExgSwiswiParamFoundType

  logical :: exg_swiswi_multi_package = .true.
  logical :: exg_swiswi_is_advanced = .false.

  character(len=16), parameter :: &
    exg_swiswi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgswiswi_inocrossstorage = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'SWISWI', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_CROSS_STORAGE', & ! tag name
    'INOCROSSSTORAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'disable cross-fluid storage Newton terms', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgswiswi_inocrossflow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'SWISWI', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_CROSS_FLOW', & ! tag name
    'INOCROSSFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'disable cross-fluid flow Newton terms', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_swiswi_param_definitions(*) = &
    [ &
    exgswiswi_inocrossstorage, &
    exgswiswi_inocrossflow &
    ]

  type(InputParamDefinitionType), parameter :: &
    exg_swiswi_aggregate_definitions(*) = &
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
    exg_swiswi_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ExgSwiswiInputModule

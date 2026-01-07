! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgGwfgweInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_gwfgwe_param_definitions
  public exg_gwfgwe_aggregate_definitions
  public exg_gwfgwe_block_definitions
  public ExgGwfgweParamFoundType
  public exg_gwfgwe_multi_package
  public exg_gwfgwe_subpackages

  type ExgGwfgweParamFoundType
  end type ExgGwfgweParamFoundType

  logical :: exg_gwfgwe_multi_package = .false.

  character(len=16), parameter :: &
    exg_gwfgwe_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exg_gwfgwe_param_definitions(*) = &
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

  type(InputParamDefinitionType), parameter :: &
    exg_gwfgwe_aggregate_definitions(*) = &
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
    exg_gwfgwe_block_definitions(*) = &
    [ &
    InputBlockDefinitionType &
    ( &
    '', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_varaible
    ) &
    ]

end module ExgGwfgweInputModule

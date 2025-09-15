! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweAdvInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_adv_param_definitions
  public gwe_adv_aggregate_definitions
  public gwe_adv_block_definitions
  public GweAdvParamFoundType
  public gwe_adv_multi_package
  public gwe_adv_subpackages

  type GweAdvParamFoundType
    logical :: scheme = .false.
    logical :: ats_percel = .false.
  end type GweAdvParamFoundType

  logical :: gwe_adv_multi_package = .false.

  character(len=16), parameter :: &
    gwe_adv_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gweadv_scheme = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'ADV', & ! subcomponent
    'OPTIONS', & ! block
    'SCHEME', & ! tag name
    'SCHEME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'advective scheme', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweadv_ats_percel = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'ADV', & ! subcomponent
    'OPTIONS', & ! block
    'ATS_PERCEL', & ! tag name
    'ATS_PERCEL', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fractional cell distance used for time step calculation', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_adv_param_definitions(*) = &
    [ &
    gweadv_scheme, &
    gweadv_ats_percel &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwe_adv_aggregate_definitions(*) = &
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
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_adv_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GweAdvInputModule

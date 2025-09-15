! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtAdvInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_adv_param_definitions
  public gwt_adv_aggregate_definitions
  public gwt_adv_block_definitions
  public GwtAdvParamFoundType
  public gwt_adv_multi_package
  public gwt_adv_subpackages

  type GwtAdvParamFoundType
    logical :: scheme = .false.
    logical :: ats_percel = .false.
  end type GwtAdvParamFoundType

  logical :: gwt_adv_multi_package = .false.

  character(len=16), parameter :: &
    gwt_adv_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtadv_scheme = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtadv_ats_percel = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwt_adv_param_definitions(*) = &
    [ &
    gwtadv_scheme, &
    gwtadv_ats_percel &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwt_adv_aggregate_definitions(*) = &
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
    gwt_adv_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwtAdvInputModule

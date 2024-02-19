! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweCndInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_cnd_param_definitions
  public gwe_cnd_aggregate_definitions
  public gwe_cnd_block_definitions
  public GweCndParamFoundType
  public gwe_cnd_multi_package

  type GweCndParamFoundType
    logical :: xt3d_off = .false.
    logical :: xt3d_rhs = .false.
    logical :: alh = .false.
    logical :: alv = .false.
    logical :: ath1 = .false.
    logical :: ath2 = .false.
    logical :: atv = .false.
    logical :: ktw = .false.
    logical :: kts = .false.
  end type GweCndParamFoundType

  logical :: gwe_cnd_multi_package = .false.

  type(InputParamDefinitionType), parameter :: &
    gwecnd_xt3d_off = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'XT3D_OFF', & ! tag name
    'XT3D_OFF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_xt3d_rhs = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'XT3D_RHS', & ! tag name
    'XT3D_RHS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_alh = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'ALH', & ! tag name
    'ALH', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_alv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'ALV', & ! tag name
    'ALV', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_ath1 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'ATH1', & ! tag name
    'ATH1', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_ath2 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'ATH2', & ! tag name
    'ATH2', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_atv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'ATV', & ! tag name
    'ATV', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_ktw = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'KTW', & ! tag name
    'KTW', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_kts = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'GRIDDATA', & ! block
    'KTS', & ! tag name
    'KTS', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_cnd_param_definitions(*) = &
    [ &
    gwecnd_xt3d_off, &
    gwecnd_xt3d_rhs, &
    gwecnd_alh, &
    gwecnd_alv, &
    gwecnd_ath1, &
    gwecnd_ath2, &
    gwecnd_atv, &
    gwecnd_ktw, &
    gwecnd_kts &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwe_cnd_aggregate_definitions(*) = &
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
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_cnd_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GweCndInputModule

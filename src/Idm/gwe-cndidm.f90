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
  public gwe_cnd_subpackages

  type GweCndParamFoundType
    logical :: xt3d_off = .false.
    logical :: xt3d_rhs = .false.
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: alh = .false.
    logical :: alv = .false.
    logical :: ath1 = .false.
    logical :: ath2 = .false.
    logical :: atv = .false.
    logical :: ktw = .false.
    logical :: kts = .false.
  end type GweCndParamFoundType

  logical :: gwe_cnd_multi_package = .false.

  character(len=16), parameter :: &
    gwe_cnd_subpackages(*) = &
    [ &
    '                ' &
    ]

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
    'deactivate xt3d', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'xt3d on right-hand side', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwecnd_export_ascii = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
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
    gwecnd_export_nc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'CND', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'longitudinal dispersivity in horizontal direction', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'longitudinal dispersivity in vertical direction', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'transverse dispersivity in horizontal direction', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'transverse dispersivity in horizontal direction', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'transverse dispersivity when flow is in vertical direction', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'thermal conductivity of the simulated fluid', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'thermal conductivity of the aquifer material', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    gwecnd_export_ascii, &
    gwecnd_export_nc, &
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

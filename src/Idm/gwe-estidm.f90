! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweEstInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_est_param_definitions
  public gwe_est_aggregate_definitions
  public gwe_est_block_definitions
  public GweEstParamFoundType
  public gwe_est_multi_package
  public gwe_est_subpackages

  type GweEstParamFoundType
    logical :: save_flows = .false.
    logical :: ord0_decay_water = .false.
    logical :: ord0_decay_solid = .false.
    logical :: rhow = .false.
    logical :: cpw = .false.
    logical :: latheatvap = .false.
    logical :: porosity = .false.
    logical :: decay_water = .false.
    logical :: decay_solid = .false.
    logical :: cps = .false.
    logical :: rhos = .false.
  end type GweEstParamFoundType

  logical :: gwe_est_multi_package = .false.

  character(len=16), parameter :: &
    gwe_est_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gweest_save_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flows to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_ord0_decay_water = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'ZERO_ORDER_DECAY_WATER', & ! tag name
    'ORD0_DECAY_WATER', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate zero-order decay in aqueous phase', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_ord0_decay_solid = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'ZERO_ORDER_DECAY_SOLID', & ! tag name
    'ORD0_DECAY_SOLID', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate zero-order decay in solid phase', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_rhow = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'DENSITY_WATER', & ! tag name
    'RHOW', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'density of water', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_cpw = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'HEAT_CAPACITY_WATER', & ! tag name
    'CPW', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'heat capacity of water', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_latheatvap = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'OPTIONS', & ! block
    'LATENT_HEAT_VAPORIZATION', & ! tag name
    'LATHEATVAP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'latent heat of vaporization', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_porosity = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'GRIDDATA', & ! block
    'POROSITY', & ! tag name
    'POROSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'porosity', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_decay_water = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY_WATER', & ! tag name
    'DECAY_WATER', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'aqueous phase decay rate coefficient', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_decay_solid = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'GRIDDATA', & ! block
    'DECAY_SOLID', & ! tag name
    'DECAY_SOLID', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'solid phase decay rate coefficient', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_cps = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'GRIDDATA', & ! block
    'HEAT_CAPACITY_SOLID', & ! tag name
    'CPS', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'heat capacity of the aquifer material', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweest_rhos = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'EST', & ! subcomponent
    'GRIDDATA', & ! block
    'DENSITY_SOLID', & ! tag name
    'RHOS', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'density of aquifer material', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_est_param_definitions(*) = &
    [ &
    gweest_save_flows, &
    gweest_ord0_decay_water, &
    gweest_ord0_decay_solid, &
    gweest_rhow, &
    gweest_cpw, &
    gweest_latheatvap, &
    gweest_porosity, &
    gweest_decay_water, &
    gweest_decay_solid, &
    gweest_cps, &
    gweest_rhos &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwe_est_aggregate_definitions(*) = &
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
    gwe_est_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GweEstInputModule

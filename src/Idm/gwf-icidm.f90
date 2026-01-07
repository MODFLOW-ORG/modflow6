! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfIcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_ic_param_definitions
  public gwf_ic_aggregate_definitions
  public gwf_ic_block_definitions
  public GwfIcParamFoundType
  public gwf_ic_multi_package
  public gwf_ic_subpackages

  type GwfIcParamFoundType
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: strt = .false.
  end type GwfIcParamFoundType

  logical :: gwf_ic_multi_package = .false.

  character(len=16), parameter :: &
    gwf_ic_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfic_export_ascii = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'IC', & ! subcomponent
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
    gwfic_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'IC', & ! subcomponent
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
    gwfic_strt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'IC', & ! subcomponent
    'GRIDDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'starting head', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_ic_param_definitions(*) = &
    [ &
    gwfic_export_ascii, &
    gwfic_export_nc, &
    gwfic_strt &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_ic_aggregate_definitions(*) = &
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
    gwf_ic_block_definitions(*) = &
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

end module GwfIcInputModule

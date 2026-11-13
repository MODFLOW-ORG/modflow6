! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlTvgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_tvg_param_definitions
  public utl_tvg_aggregate_definitions
  public utl_tvg_block_definitions
  public UtlTvgParamFoundType
  public utl_tvg_multi_package
  public utl_tvg_is_advanced
  public utl_tvg_subpackages

  type UtlTvgParamFoundType
    logical :: readarraygrid = .false.
    logical :: auxiliary = .false.
    logical :: iprpak = .false.
    logical :: export_nc = .false.
    logical :: auxvar = .false.
  end type UtlTvgParamFoundType

  logical :: utl_tvg_multi_package = .true.
  logical :: utl_tvg_is_advanced = .false.

  character(len=16), parameter :: &
    utl_tvg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltvg_readarraygrid = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVG', & ! subcomponent
    'OPTIONS', & ! block
    'READARRAYGRID', & ! tag name
    'READARRAYGRID', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use grid-array input', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvg_auxiliary = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVG', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvg_iprpak = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVG', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvg_export_nc = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVG', & ! subcomponent
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
    utltvg_auxvar = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvg_param_definitions(*) = &
    [ &
    utltvg_readarraygrid, &
    utltvg_auxiliary, &
    utltvg_iprpak, &
    utltvg_export_nc, &
    utltvg_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    utl_tvg_aggregate_definitions(*) = &
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
    utl_tvg_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
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

end module UtlTvgInputModule

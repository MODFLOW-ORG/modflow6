! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlTvkgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_tvkg_param_definitions
  public utl_tvkg_aggregate_definitions
  public utl_tvkg_block_definitions
  public UtlTvkgParamFoundType
  public utl_tvkg_multi_package
  public utl_tvkg_is_advanced
  public utl_tvkg_subpackages

  type UtlTvkgParamFoundType
    logical :: readarraygrid = .false.
    logical :: print_input = .false.
    logical :: tas_filerecord = .false.
    logical :: tas6 = .false.
    logical :: filein = .false.
    logical :: tas6_filename = .false.
    logical :: export_nc = .false.
    logical :: k = .false.
    logical :: k22 = .false.
    logical :: k33 = .false.
  end type UtlTvkgParamFoundType

  logical :: utl_tvkg_multi_package = .false.
  logical :: utl_tvkg_is_advanced = .false.

  character(len=16), parameter :: &
    utl_tvkg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltvkg_readarraygrid = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'READARRAYGRID', & ! tag name
    'READARRAYGRID', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use array-based grid input', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
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
    utltvkg_tas_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'TAS_FILERECORD', & ! tag name
    'TAS_FILERECORD', & ! fortran variable
    'RECORD TAS6 FILEIN TAS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_tas6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6', & ! tag name
    'TAS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'time-array series keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_tas6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6_FILENAME', & ! tag name
    'TAS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_export_nc = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
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
    utltvkg_k = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'PERIOD', & ! block
    'K', & ! tag name
    'K', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity (L/T)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_k22 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'PERIOD', & ! block
    'K22', & ! tag name
    'K22', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity of second ellipsoid axis (L/T)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvkg_k33 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVKG', & ! subcomponent
    'PERIOD', & ! block
    'K33', & ! tag name
    'K33', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'hydraulic conductivity of third ellipsoid axis (L/T)', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvkg_param_definitions(*) = &
    [ &
    utltvkg_readarraygrid, &
    utltvkg_print_input, &
    utltvkg_tas_filerecord, &
    utltvkg_tas6, &
    utltvkg_filein, &
    utltvkg_tas6_filename, &
    utltvkg_export_nc, &
    utltvkg_k, &
    utltvkg_k22, &
    utltvkg_k33 &
    ]

  type(InputParamDefinitionType), parameter :: &
    utl_tvkg_aggregate_definitions(*) = &
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
    utl_tvkg_block_definitions(*) = &
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

end module UtlTvkgInputModule

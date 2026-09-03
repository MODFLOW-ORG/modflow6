! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlTvkInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_tvk_param_definitions
  public utl_tvk_aggregate_definitions
  public utl_tvk_block_definitions
  public UtlTvkParamFoundType
  public utl_tvk_multi_package
  public utl_tvk_subpackages

  type UtlTvkParamFoundType
    logical :: print_input = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: cellid = .false.
    logical :: k_in = .false.
    logical :: k22_in = .false.
    logical :: k33_in = .false.
  end type UtlTvkParamFoundType

  logical :: utl_tvk_multi_package = .false.

  character(len=16), parameter :: &
    utl_tvk_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltvk_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
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
    utltvk_ts_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
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
    utltvk_ts6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvk_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
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
    utltvk_ts6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
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
    utltvk_cellid = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvk_k_in = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'K', & ! tag name
    'K_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'hydraulic conductivity (L/T)', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvk_k22_in = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'K22', & ! tag name
    'K22_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'hydraulic conductivity of second ellipsoid axis (L/T)', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvk_k33_in = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'K33', & ! tag name
    'K33_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'hydraulic conductivity of third ellipsoid axis (L/T)', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvk_param_definitions(*) = &
    [ &
    utltvk_print_input, &
    utltvk_ts_filerecord, &
    utltvk_ts6, &
    utltvk_filein, &
    utltvk_ts6_filename, &
    utltvk_cellid, &
    utltvk_k_in, &
    utltvk_k22_in, &
    utltvk_k33_in &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltvk_perioddata = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY CELLID TVKSETTING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvk_tvksetting = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVK', & ! subcomponent
    'PERIOD', & ! block
    'TVKSETTING', & ! tag name
    'TVKSETTING', & ! fortran variable
    'KEYSTRING K K22 K33', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvk_aggregate_definitions(*) = &
    [ &
    utltvk_perioddata, &
    utltvk_tvksetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_tvk_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module UtlTvkInputModule

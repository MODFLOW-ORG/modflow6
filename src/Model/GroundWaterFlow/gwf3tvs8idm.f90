! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlTvsInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_tvs_param_definitions
  public utl_tvs_aggregate_definitions
  public utl_tvs_block_definitions
  public UtlTvsParamFoundType
  public utl_tvs_multi_package

  type UtlTvsParamFoundType
    logical :: disablestochg = .false.
    logical :: print_input = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: cellid = .false.
    logical :: tvssetting = .false.
    logical :: ss = .false.
    logical :: sy = .false.
  end type UtlTvsParamFoundType

  logical :: utl_tvs_multi_package = .false.

  type(InputParamDefinitionType), parameter :: &
    utltvs_disablestochg = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'DISABLE_STORAGE_CHANGE_INTEGRATION', & ! tag name
    'DISABLESTOCHG', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_ts_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_ts6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_ts6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_cellid = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_tvssetting = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'PERIOD', & ! block
    'TVSSETTING', & ! tag name
    'TVSSETTING', & ! fortran variable
    'KEYSTRING SS SY', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_ss = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'PERIOD', & ! block
    'SS', & ! tag name
    'SS', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltvs_sy = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'PERIOD', & ! block
    'SY', & ! tag name
    'SY', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvs_param_definitions(*) = &
    [ &
    utltvs_disablestochg, &
    utltvs_print_input, &
    utltvs_ts_filerecord, &
    utltvs_ts6, &
    utltvs_filein, &
    utltvs_ts6_filename, &
    utltvs_cellid, &
    utltvs_tvssetting, &
    utltvs_ss, &
    utltvs_sy &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltvs_perioddata = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVS', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY CELLID TVSSETTING', & ! type
    '', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tvs_aggregate_definitions(*) = &
    [ &
    utltvs_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_tvs_block_definitions(*) = &
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

end module UtlTvsInputModule

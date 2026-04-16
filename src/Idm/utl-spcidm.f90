! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlSpcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_spc_param_definitions
  public utl_spc_aggregate_definitions
  public utl_spc_block_definitions
  public UtlSpcParamFoundType
  public utl_spc_multi_package
  public utl_spc_subpackages

  type UtlSpcParamFoundType
    logical :: print_input = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: maxbound = .false.
    logical :: bndno = .false.
    logical :: concentration = .false.
    logical :: temperature = .false.
  end type UtlSpcParamFoundType

  logical :: utl_spc_multi_package = .true.

  character(len=16), parameter :: &
    utl_spc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlspc_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
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
    utlspc_ts_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
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
    utlspc_ts6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
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
    utlspc_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
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
    utlspc_ts6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
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
    utlspc_maxbound = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of ssm entries', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlspc_bndno = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'PERIOD', & ! block
    'BNDNO', & ! tag name
    'BNDNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'bound number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlspc_concentration = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'boundary concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlspc_temperature = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'boundary temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_spc_param_definitions(*) = &
    [ &
    utlspc_print_input, &
    utlspc_ts_filerecord, &
    utlspc_ts6, &
    utlspc_filein, &
    utlspc_ts6_filename, &
    utlspc_maxbound, &
    utlspc_bndno, &
    utlspc_concentration, &
    utlspc_temperature &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlspc_perioddata = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY BNDNO SPCSETTING', & ! type
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
    utlspc_spcsetting = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPC', & ! subcomponent
    'PERIOD', & ! block
    'SPCSETTING', & ! tag name
    'SPCSETTING', & ! fortran variable
    'KEYSTRING CONCENTRATION TEMPERATURE', & ! type
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
    utl_spc_aggregate_definitions(*) = &
    [ &
    utlspc_perioddata, &
    utlspc_spcsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_spc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
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

end module UtlSpcInputModule

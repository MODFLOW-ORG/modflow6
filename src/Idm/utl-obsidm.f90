! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlObsInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_obs_param_definitions
  public utl_obs_aggregate_definitions
  public utl_obs_block_definitions
  public UtlObsParamFoundType
  public utl_obs_multi_package
  public utl_obs_subpackages

  type UtlObsParamFoundType
    logical :: digits = .false.
    logical :: print_input = .false.
    logical :: fileout = .false.
    logical :: obs_output_fname = .false.
    logical :: binary = .false.
    logical :: obsname = .false.
    logical :: obstype = .false.
    logical :: id = .false.
    logical :: id2 = .false.
  end type UtlObsParamFoundType

  logical :: utl_obs_multi_package = .true.

  character(len=16), parameter :: &
    utl_obs_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlobs_digits = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'OPTIONS', & ! block
    'DIGITS', & ! tag name
    'DIGITS', & ! fortran variable
    'INTEGER', & ! type
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
    utlobs_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
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
    utlobs_fileout = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
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
    utlobs_obs_output_fname = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'OBS_OUTPUT_FILE_NAME', & ! tag name
    'OBS_OUTPUT_FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlobs_binary = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'BINARY', & ! tag name
    'BINARY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlobs_obsname = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'OBSNAME', & ! tag name
    'OBSNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'observation name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlobs_obstype = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'OBSTYPE', & ! tag name
    'OBSTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'observation type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlobs_id = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'ID', & ! tag name
    'ID', & ! fortran variable
    'STRING', & ! type
    'LINELENGTH', & ! shape
    'time', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlobs_id2 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'ID2', & ! tag name
    'ID2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'time', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_obs_param_definitions(*) = &
    [ &
    utlobs_digits, &
    utlobs_print_input, &
    utlobs_fileout, &
    utlobs_obs_output_fname, &
    utlobs_binary, &
    utlobs_obsname, &
    utlobs_obstype, &
    utlobs_id, &
    utlobs_id2 &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlobs_continuous = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'OBS', & ! subcomponent
    'CONTINUOUS', & ! block
    'CONTINUOUS', & ! tag name
    'CONTINUOUS', & ! fortran variable
    'RECARRAY OBSNAME OBSTYPE ID ID2', & ! type
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
    utl_obs_aggregate_definitions(*) = &
    [ &
    utlobs_continuous &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_obs_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CONTINUOUS', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module UtlObsInputModule

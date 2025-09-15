! ** Do Not Modify! MODFLOW 6 system generated file. **
module SimTdisInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public sim_tdis_param_definitions
  public sim_tdis_aggregate_definitions
  public sim_tdis_block_definitions
  public SimTdisParamFoundType
  public sim_tdis_multi_package
  public sim_tdis_subpackages

  type SimTdisParamFoundType
    logical :: time_units = .false.
    logical :: start_date_time = .false.
    logical :: ats_filerecord = .false.
    logical :: ats6 = .false.
    logical :: filein = .false.
    logical :: ats6_filename = .false.
    logical :: nper = .false.
    logical :: perlen = .false.
    logical :: nstp = .false.
    logical :: tsmult = .false.
  end type SimTdisParamFoundType

  logical :: sim_tdis_multi_package = .false.

  character(len=16), parameter :: &
    sim_tdis_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    simtdis_time_units = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_UNITS', & ! tag name
    'TIME_UNITS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'time unit', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_start_date_time = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'START_DATE_TIME', & ! tag name
    'START_DATE_TIME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'starting date and time', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_ats_filerecord = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'ATS_FILERECORD', & ! tag name
    'ATS_FILERECORD', & ! fortran variable
    'RECORD ATS6 FILEIN ATS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_ats6 = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'ATS6', & ! tag name
    'ATS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'ats keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_filein = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_ats6_filename = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'OPTIONS', & ! block
    'ATS6_FILENAME', & ! tag name
    'ATS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of adaptive time series information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_nper = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NPER', & ! tag name
    'NPER', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of stress periods', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_perlen = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'PERIODDATA', & ! block
    'PERLEN', & ! tag name
    'PERLEN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'length of stress period', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_nstp = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'PERIODDATA', & ! block
    'NSTP', & ! tag name
    'NSTP', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of time steps', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    simtdis_tsmult = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'PERIODDATA', & ! block
    'TSMULT', & ! tag name
    'TSMULT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'number of time steps', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    sim_tdis_param_definitions(*) = &
    [ &
    simtdis_time_units, &
    simtdis_start_date_time, &
    simtdis_ats_filerecord, &
    simtdis_ats6, &
    simtdis_filein, &
    simtdis_ats6_filename, &
    simtdis_nper, &
    simtdis_perlen, &
    simtdis_nstp, &
    simtdis_tsmult &
    ]

  type(InputParamDefinitionType), parameter :: &
    simtdis_perioddata = InputParamDefinitionType &
    ( &
    'SIM', & ! component
    'TDIS', & ! subcomponent
    'PERIODDATA', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY PERLEN NSTP TSMULT', & ! type
    '', & ! shape
    'stress period time information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    sim_tdis_aggregate_definitions(*) = &
    [ &
    simtdis_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    sim_tdis_block_definitions(*) = &
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
    'PERIODDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module SimTdisInputModule

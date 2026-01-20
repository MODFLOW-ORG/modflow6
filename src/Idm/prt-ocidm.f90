! ** Do Not Modify! MODFLOW 6 system generated file. **
module PrtOcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public prt_oc_param_definitions
  public prt_oc_aggregate_definitions
  public prt_oc_block_definitions
  public PrtOcParamFoundType
  public prt_oc_multi_package
  public prt_oc_subpackages

  type PrtOcParamFoundType
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: trackfilerec = .false.
    logical :: track = .false.
    logical :: trackfile = .false.
    logical :: trackcsvfilerec = .false.
    logical :: trackcsv = .false.
    logical :: trackcsvfile = .false.
    logical :: track_release = .false.
    logical :: track_exit = .false.
    logical :: track_subf_exit = .false.
    logical :: track_timestep = .false.
    logical :: track_terminate = .false.
    logical :: track_weaksink = .false.
    logical :: track_usertime = .false.
    logical :: track_dropped = .false.
    logical :: ttimesrec = .false.
    logical :: track_times = .false.
    logical :: times = .false.
    logical :: ttimesfilerec = .false.
    logical :: track_timesfile = .false.
    logical :: timesfile = .false.
    logical :: dev_dump_evtrace = .false.
    logical :: ntracktimes = .false.
    logical :: time = .false.
    logical :: saverecord = .false.
    logical :: save = .false.
    logical :: printrecord = .false.
    logical :: print = .false.
    logical :: rtype = .false.
    logical :: ocsetting = .false.
    logical :: all = .false.
    logical :: first = .false.
    logical :: last = .false.
    logical :: frequency = .false.
    logical :: steps = .false.
  end type PrtOcParamFoundType

  logical :: prt_oc_multi_package = .false.

  character(len=16), parameter :: &
    prt_oc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtoc_budfilerec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET_FILERECORD', & ! tag name
    'BUDFILEREC', & ! fortran variable
    'RECORD BUDGET FILEOUT BUDGETFILE', & ! type
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
    prtoc_budget = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET', & ! tag name
    'BUDGET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_fileout = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
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
    prtoc_budgetfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETFILE', & ! tag name
    'BUDGETFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_budcsvfilerec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV_FILERECORD', & ! tag name
    'BUDCSVFILEREC', & ! fortran variable
    'RECORD BUDGETCSV FILEOUT BUDGETCSVFILE', & ! type
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
    prtoc_budgetcsv = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV', & ! tag name
    'BUDGETCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_budgetcsvfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSVFILE', & ! tag name
    'BUDGETCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_trackfilerec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_FILERECORD', & ! tag name
    'TRACKFILEREC', & ! fortran variable
    'RECORD TRACK FILEOUT TRACKFILE', & ! type
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
    prtoc_track = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK', & ! tag name
    'TRACK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_trackfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACKFILE', & ! tag name
    'TRACKFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_trackcsvfilerec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACKCSV_FILERECORD', & ! tag name
    'TRACKCSVFILEREC', & ! fortran variable
    'RECORD TRACKCSV FILEOUT TRACKCSVFILE', & ! type
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
    prtoc_trackcsv = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACKCSV', & ! tag name
    'TRACKCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_trackcsvfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACKCSVFILE', & ! tag name
    'TRACKCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_release = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_RELEASE', & ! tag name
    'TRACK_RELEASE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track release', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_exit = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_EXIT', & ! tag name
    'TRACK_EXIT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track domain exits', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_subf_exit = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_SUBFEATURE_EXIT', & ! tag name
    'TRACK_SUBF_EXIT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track cell exits', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_timestep = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TIMESTEP', & ! tag name
    'TRACK_TIMESTEP', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track timestep ends', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_terminate = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TERMINATE', & ! tag name
    'TRACK_TERMINATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track termination', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_weaksink = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_WEAKSINK', & ! tag name
    'TRACK_WEAKSINK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track weaksink exits', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_usertime = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_USERTIME', & ! tag name
    'TRACK_USERTIME', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track user-specified times', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_track_dropped = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_DROPPED', & ! tag name
    'TRACK_DROPPED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'track drops to water table', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_ttimesrec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TIMESRECORD', & ! tag name
    'TTIMESREC', & ! fortran variable
    'RECORD TRACK_TIMES TIMES', & ! type
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
    prtoc_track_times = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TIMES', & ! tag name
    'TRACK_TIMES', & ! fortran variable
    'KEYWORD', & ! type
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
    prtoc_times = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TIMES', & ! tag name
    'TIMES', & ! fortran variable
    'DOUBLE1D', & ! type
    'ANY1D', & ! shape
    'tracking times', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_ttimesfilerec = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TIMESFILERECORD', & ! tag name
    'TTIMESFILEREC', & ! fortran variable
    'RECORD TRACK_TIMESFILE TIMESFILE', & ! type
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
    prtoc_track_timesfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_TIMESFILE', & ! tag name
    'TRACK_TIMESFILE', & ! fortran variable
    'KEYWORD', & ! type
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
    prtoc_timesfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'TIMESFILE', & ! tag name
    'TIMESFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_dev_dump_evtrace = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_DUMP_EVENT_TRACE', & ! tag name
    'DEV_DUMP_EVTRACE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print particle tracking events', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_ntracktimes = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'DIMENSIONS', & ! block
    'NTRACKTIMES', & ! tag name
    'NTRACKTIMES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of particle tracking times', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_time = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'TRACKTIMES', & ! block
    'TIME', & ! tag name
    'TIME', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'release time', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_saverecord = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'SAVERECORD', & ! tag name
    'SAVERECORD', & ! fortran variable
    'RECORD SAVE RTYPE OCSETTING', & ! type
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
    prtoc_save = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'SAVE', & ! tag name
    'SAVE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_printrecord = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'PRINTRECORD', & ! tag name
    'PRINTRECORD', & ! fortran variable
    'RECORD PRINT RTYPE OCSETTING', & ! type
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
    prtoc_print = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'PRINT', & ! tag name
    'PRINT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_rtype = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'RTYPE', & ! tag name
    'RTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'record type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtoc_ocsetting = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'OCSETTING', & ! tag name
    'OCSETTING', & ! fortran variable
    'KEYSTRING ALL FIRST LAST FREQUENCY STEPS', & ! type
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
    prtoc_all = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'ALL', & ! tag name
    'ALL', & ! fortran variable
    'KEYWORD', & ! type
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
    prtoc_first = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'FIRST', & ! tag name
    'FIRST', & ! fortran variable
    'KEYWORD', & ! type
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
    prtoc_last = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'LAST', & ! tag name
    'LAST', & ! fortran variable
    'KEYWORD', & ! type
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
    prtoc_frequency = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'FREQUENCY', & ! tag name
    'FREQUENCY', & ! fortran variable
    'INTEGER', & ! type
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
    prtoc_steps = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'STEPS', & ! tag name
    'STEPS', & ! fortran variable
    'INTEGER1D', & ! type
    '<NSTP', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prt_oc_param_definitions(*) = &
    [ &
    prtoc_budfilerec, &
    prtoc_budget, &
    prtoc_fileout, &
    prtoc_budgetfile, &
    prtoc_budcsvfilerec, &
    prtoc_budgetcsv, &
    prtoc_budgetcsvfile, &
    prtoc_trackfilerec, &
    prtoc_track, &
    prtoc_trackfile, &
    prtoc_trackcsvfilerec, &
    prtoc_trackcsv, &
    prtoc_trackcsvfile, &
    prtoc_track_release, &
    prtoc_track_exit, &
    prtoc_track_subf_exit, &
    prtoc_track_timestep, &
    prtoc_track_terminate, &
    prtoc_track_weaksink, &
    prtoc_track_usertime, &
    prtoc_track_dropped, &
    prtoc_ttimesrec, &
    prtoc_track_times, &
    prtoc_times, &
    prtoc_ttimesfilerec, &
    prtoc_track_timesfile, &
    prtoc_timesfile, &
    prtoc_dev_dump_evtrace, &
    prtoc_ntracktimes, &
    prtoc_time, &
    prtoc_saverecord, &
    prtoc_save, &
    prtoc_printrecord, &
    prtoc_print, &
    prtoc_rtype, &
    prtoc_ocsetting, &
    prtoc_all, &
    prtoc_first, &
    prtoc_last, &
    prtoc_frequency, &
    prtoc_steps &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtoc_tracktimes = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'OC', & ! subcomponent
    'TRACKTIMES', & ! block
    'TRACKTIMES', & ! tag name
    'TRACKTIMES', & ! fortran variable
    'RECARRAY TIME', & ! type
    'NTRACKTIMES', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prt_oc_aggregate_definitions(*) = &
    [ &
    prtoc_tracktimes &
    ]

  type(InputBlockDefinitionType), parameter :: &
    prt_oc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'TRACKTIMES', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module PrtOcInputModule

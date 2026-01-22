! ** Do Not Modify! MODFLOW 6 system generated file. **
module PrtPrpInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public prt_prp_param_definitions
  public prt_prp_aggregate_definitions
  public prt_prp_block_definitions
  public PrtPrpParamFoundType
  public prt_prp_multi_package
  public prt_prp_subpackages

  type PrtPrpParamFoundType
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iexmeth = .false.
    logical :: extol = .false.
    logical :: localz = .false.
    logical :: extend = .false.
    logical :: track_filerecord = .false.
    logical :: track = .false.
    logical :: fileout = .false.
    logical :: trackfile = .false.
    logical :: trackcsvfr = .false.
    logical :: trackcsv = .false.
    logical :: trackcsvfile = .false.
    logical :: stoptime = .false.
    logical :: stoptraveltime = .false.
    logical :: istopweaksink = .false.
    logical :: istopzone = .false.
    logical :: drape = .false.
    logical :: releasetr = .false.
    logical :: release_times = .false.
    logical :: times = .false.
    logical :: release_timesfr = .false.
    logical :: release_timesfn = .false.
    logical :: timesfile = .false.
    logical :: idrymeth = .false.
    logical :: frctrn = .false.
    logical :: rttol = .false.
    logical :: rtfreq = .false.
    logical :: ichkmeth = .false.
    logical :: icycwin = .false.
    logical :: nreleasepts = .false.
    logical :: nreleasetimes = .false.
    logical :: irptno = .false.
    logical :: cellid = .false.
    logical :: xrpt = .false.
    logical :: yrpt = .false.
    logical :: zrpt = .false.
    logical :: boundname = .false.
    logical :: time = .false.
    logical :: releasesetting = .false.
    logical :: all = .false.
    logical :: first = .false.
    logical :: last = .false.
    logical :: frequency = .false.
    logical :: steps = .false.
    logical :: fraction = .false.
  end type PrtPrpParamFoundType

  logical :: prt_prp_multi_package = .true.

  character(len=16), parameter :: &
    prt_prp_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtprp_boundnames = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
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
    prtprp_iprpak = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_iexmeth = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_EXIT_SOLVE_METHOD', & ! tag name
    'IEXMETH', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'exit solve method', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_extol = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'EXIT_SOLVE_TOLERANCE', & ! tag name
    'EXTOL', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'exit solve tolerance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_localz = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'LOCAL_Z', & ! tag name
    'LOCALZ', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'whether to use local z coordinates', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_extend = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'EXTEND_TRACKING', & ! tag name
    'EXTEND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'whether to extend tracking beyond the end of the simulation', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_track_filerecord = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'TRACK_FILERECORD', & ! tag name
    'TRACK_FILERECORD', & ! fortran variable
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
    prtprp_track = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_fileout = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_trackfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_trackcsvfr = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'TRACKCSV_FILERECORD', & ! tag name
    'TRACKCSVFR', & ! fortran variable
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
    prtprp_trackcsv = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_trackcsvfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_stoptime = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'STOPTIME', & ! tag name
    'STOPTIME', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'stop time', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_stoptraveltime = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'STOPTRAVELTIME', & ! tag name
    'STOPTRAVELTIME', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'stop travel time', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_istopweaksink = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'STOP_AT_WEAK_SINK', & ! tag name
    'ISTOPWEAKSINK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'stop at weak sink', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_istopzone = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'ISTOPZONE', & ! tag name
    'ISTOPZONE', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'stop zone number', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_drape = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'DRAPE', & ! tag name
    'DRAPE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'drape', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_releasetr = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIMESRECORD', & ! tag name
    'RELEASETR', & ! fortran variable
    'RECORD RELEASE_TIMES TIMES', & ! type
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
    prtprp_release_times = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIMES', & ! tag name
    'RELEASE_TIMES', & ! fortran variable
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
    prtprp_times = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'TIMES', & ! tag name
    'TIMES', & ! fortran variable
    'DOUBLE1D', & ! type
    'ANY1D', & ! shape
    'release times', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_release_timesfr = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIMESFILERECORD', & ! tag name
    'RELEASE_TIMESFR', & ! fortran variable
    'RECORD RELEASE_TIMESFILE TIMESFILE', & ! type
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
    prtprp_release_timesfn = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIMESFILE', & ! tag name
    'RELEASE_TIMESFN', & ! fortran variable
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
    prtprp_timesfile = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_idrymeth = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'DRY_TRACKING_METHOD', & ! tag name
    'IDRYMETH', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'what to do in dry-but-active cells', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_frctrn = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_FORCETERNARY', & ! tag name
    'FRCTRN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'force ternary tracking method', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_rttol = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIME_TOLERANCE', & ! tag name
    'RTTOL', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'release time coincidence tolerance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_rtfreq = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'RELEASE_TIME_FREQUENCY', & ! tag name
    'RTFREQ', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'release time frequency', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_ichkmeth = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'COORDINATE_CHECK_METHOD', & ! tag name
    'ICHKMETH', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'coordinate checking method', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_icycwin = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_CYCLE_DETECTION_WINDOW', & ! tag name
    'ICYCWIN', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cycle detection window size', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_nreleasepts = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'DIMENSIONS', & ! block
    'NRELEASEPTS', & ! tag name
    'NRELEASEPTS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of particle release points', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_nreleasetimes = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'DIMENSIONS', & ! block
    'NRELEASETIMES', & ! tag name
    'NRELEASETIMES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of particle release times', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_irptno = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IRPTNO', & ! tag name
    'IRPTNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'PRP id number for release point', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_cellid = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
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
    prtprp_xrpt = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'XRPT', & ! tag name
    'XRPT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x coordinate of release point', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_yrpt = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'YRPT', & ! tag name
    'YRPT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y coordinate of release point', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_zrpt = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'ZRPT', & ! tag name
    'ZRPT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'z coordinate of release point', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_boundname = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'release point name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_time = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'RELEASETIMES', & ! block
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
    prtprp_releasesetting = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PERIOD', & ! block
    'RELEASESETTING', & ! tag name
    'RELEASESETTING', & ! fortran variable
    'KEYSTRING ALL FIRST LAST FREQUENCY STEPS FRACTION', & ! type
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
    prtprp_all = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_first = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_last = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_frequency = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_steps = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
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
    prtprp_fraction = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PERIOD', & ! block
    'FRACTION', & ! tag name
    'FRACTION', & ! fortran variable
    'DOUBLE1D', & ! type
    '<NSTP', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prt_prp_param_definitions(*) = &
    [ &
    prtprp_boundnames, &
    prtprp_iprpak, &
    prtprp_iexmeth, &
    prtprp_extol, &
    prtprp_localz, &
    prtprp_extend, &
    prtprp_track_filerecord, &
    prtprp_track, &
    prtprp_fileout, &
    prtprp_trackfile, &
    prtprp_trackcsvfr, &
    prtprp_trackcsv, &
    prtprp_trackcsvfile, &
    prtprp_stoptime, &
    prtprp_stoptraveltime, &
    prtprp_istopweaksink, &
    prtprp_istopzone, &
    prtprp_drape, &
    prtprp_releasetr, &
    prtprp_release_times, &
    prtprp_times, &
    prtprp_release_timesfr, &
    prtprp_release_timesfn, &
    prtprp_timesfile, &
    prtprp_idrymeth, &
    prtprp_frctrn, &
    prtprp_rttol, &
    prtprp_rtfreq, &
    prtprp_ichkmeth, &
    prtprp_icycwin, &
    prtprp_nreleasepts, &
    prtprp_nreleasetimes, &
    prtprp_irptno, &
    prtprp_cellid, &
    prtprp_xrpt, &
    prtprp_yrpt, &
    prtprp_zrpt, &
    prtprp_boundname, &
    prtprp_time, &
    prtprp_releasesetting, &
    prtprp_all, &
    prtprp_first, &
    prtprp_last, &
    prtprp_frequency, &
    prtprp_steps, &
    prtprp_fraction &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtprp_packagedata = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IRPTNO CELLID XRPT YRPT ZRPT BOUNDNAME', & ! type
    'NRELEASEPTS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_releasetimes = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'RELEASETIMES', & ! block
    'RELEASETIMES', & ! tag name
    'RELEASETIMES', & ! fortran variable
    'RECARRAY TIME', & ! type
    'NRELEASETIMES', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtprp_perioddata = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'PRP', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY RELEASESETTING', & ! type
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
    prt_prp_aggregate_definitions(*) = &
    [ &
    prtprp_packagedata, &
    prtprp_releasetimes, &
    prtprp_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    prt_prp_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
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
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'RELEASETIMES', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module PrtPrpInputModule

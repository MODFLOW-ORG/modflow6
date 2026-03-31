! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfMawInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_maw_param_definitions
  public gwf_maw_aggregate_definitions
  public gwf_maw_block_definitions
  public GwfMawParamFoundType
  public gwf_maw_multi_package
  public gwf_maw_subpackages

  type GwfMawParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: print_head = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: head_filerecord = .false.
    logical :: head = .false.
    logical :: headfile = .false.
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: no_well_storage = .false.
    logical :: flow_correction = .false.
    logical :: inonvert = .false.
    logical :: flowing_wells = .false.
    logical :: shutdown_theta = .false.
    logical :: shutdown_kappa = .false.
    logical :: mfrcsvfilerec = .false.
    logical :: maw_flw_red_csv = .false.
    logical :: mfrcsvfile = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: ieffradopt = .false.
    logical :: nmawwells = .false.
    logical :: packagedata_ifno = .false.
    logical :: radius = .false.
    logical :: bottom = .false.
    logical :: strt = .false.
    logical :: condeqn = .false.
    logical :: ngwfnodes = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: conndata_ifno = .false.
    logical :: icon = .false.
    logical :: cellid = .false.
    logical :: scrn_top = .false.
    logical :: scrn_bot = .false.
    logical :: hk_skin = .false.
    logical :: radius_skin = .false.
    logical :: angledata_ifno = .false.
    logical :: angledata_icon = .false.
    logical :: angle = .false.
    logical :: connlen = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: flwwellrecord = .false.
    logical :: flowing_well = .false.
    logical :: fwelev = .false.
    logical :: fwcond = .false.
    logical :: fwrlen = .false.
    logical :: rate = .false.
    logical :: well_head = .false.
    logical :: head_limit = .false.
    logical :: shutoffrecord = .false.
    logical :: shut_off = .false.
    logical :: minrate = .false.
    logical :: maxrate = .false.
    logical :: ratescalingrec = .false.
    logical :: rate_scaling = .false.
    logical :: pump_elevation = .false.
    logical :: scaling_length = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwfMawParamFoundType

  logical :: gwf_maw_multi_package = .true.

  character(len=16), parameter :: &
    gwf_maw_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_print_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_HEAD', & ! tag name
    'PRINT_HEAD', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated heads to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save well flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_head_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'HEAD_FILERECORD', & ! tag name
    'HEAD_FILERECORD', & ! fortran variable
    'RECORD HEAD FILEOUT HEADFILE', & ! type
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
    gwfmaw_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'HEAD', & ! tag name
    'HEAD', & ! fortran variable
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
    gwfmaw_headfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'HEADFILE', & ! tag name
    'HEADFILE', & ! fortran variable
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
    gwfmaw_budfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_no_well_storage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'NO_WELL_STORAGE', & ! tag name
    'NO_WELL_STORAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'deactivate well storage', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_flow_correction = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_CORRECTION', & ! tag name
    'FLOW_CORRECTION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate flow correction', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_inonvert = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'NON_VERTICAL_WELLS', & ! tag name
    'INONVERT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate non-vertical well connections', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_flowing_wells = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'FLOWING_WELLS', & ! tag name
    'FLOWING_WELLS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate flowing wells', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_shutdown_theta = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'SHUTDOWN_THETA', & ! tag name
    'SHUTDOWN_THETA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'shutdown theta', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_shutdown_kappa = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'SHUTDOWN_KAPPA', & ! tag name
    'SHUTDOWN_KAPPA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'shutdown kappa', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_mfrcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'MFRCSV_FILERECORD', & ! tag name
    'MFRCSVFILEREC', & ! fortran variable
    'RECORD MAW_FLOW_REDUCE_CSV FILEOUT MFRCSVFILE', & ! type
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
    gwfmaw_maw_flw_red_csv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'MAW_FLOW_REDUCE_CSV', & ! tag name
    'MAW_FLW_RED_CSV', & ! fortran variable
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
    gwfmaw_mfrcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'MFRCSVFILE', & ! tag name
    'MFRCSVFILE', & ! fortran variable
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
    gwfmaw_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
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
    gwfmaw_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
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
    gwfmaw_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'MOVER', & ! tag name
    'MOVER', & ! fortran variable
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
    gwfmaw_ieffradopt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_PEACEMAN_EFFECTIVE_RADIUS', & ! tag name
    'IEFFRADOPT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use Peaceman effective radius for structured grids', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_nmawwells = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'DIMENSIONS', & ! block
    'NMAWWELLS', & ! tag name
    'NMAWWELLS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of MAW wells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_radius = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RADIUS', & ! tag name
    'RADIUS', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well radius', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_bottom = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOTTOM', & ! tag name
    'BOTTOM', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well bottom', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_strt = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting head', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_condeqn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'CONDEQN', & ! tag name
    'CONDEQN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'conductance equation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_ngwfnodes = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NGWFNODES', & ! tag name
    'NGWFNODES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of connected GWF cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_aux = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'AUX', & ! tag name
    'AUX', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_conndata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IFNO', & ! tag name
    'CONNDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_icon = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'ICON', & ! tag name
    'ICON', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'connection number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
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
    gwfmaw_scrn_top = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'SCRN_TOP', & ! tag name
    'SCRN_TOP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'screen top', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_scrn_bot = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'SCRN_BOT', & ! tag name
    'SCRN_BOT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'screen bottom', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_hk_skin = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'HK_SKIN', & ! tag name
    'HK_SKIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'skin data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_radius_skin = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'RADIUS_SKIN', & ! tag name
    'RADIUS_SKIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'skin radius', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_angledata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'ANGLEDATA', & ! block
    'IFNO', & ! tag name
    'ANGLEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_angledata_icon = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'ANGLEDATA', & ! block
    'ICON', & ! tag name
    'ANGLEDATA_ICON', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'connection number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_angle = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'ANGLEDATA', & ! block
    'ANGLE', & ! tag name
    'ANGLE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'tilt angle from vertical', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_connlen = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'ANGLEDATA', & ! block
    'CONN_LENGTH', & ! tag name
    'CONNLEN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'in-cell screen length', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_status = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_flwwellrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'FLOWING_WELLRECORD', & ! tag name
    'FLWWELLRECORD', & ! fortran variable
    'RECORD FLOWING_WELL FWELEV FWCOND FWRLEN', & ! type
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
    gwfmaw_flowing_well = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'FLOWING_WELL', & ! tag name
    'FLOWING_WELL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'well is a flowing well', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_fwelev = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'FWELEV', & ! tag name
    'FWELEV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'flowing well elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_fwcond = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'FWCOND', & ! tag name
    'FWCOND', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well flowing well conductance', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_fwrlen = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'FWRLEN', & ! tag name
    'FWRLEN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'flowing well reduction length', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_rate = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well pumping rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_well_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'WELL_HEAD', & ! tag name
    'WELL_HEAD', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well head', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_head_limit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'HEAD_LIMIT', & ! tag name
    'HEAD_LIMIT', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'head limit', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_shutoffrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'SHUTOFFRECORD', & ! tag name
    'SHUTOFFRECORD', & ! fortran variable
    'RECORD SHUT_OFF MINRATE MAXRATE', & ! type
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
    gwfmaw_shut_off = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'SHUT_OFF', & ! tag name
    'SHUT_OFF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'shut off well', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_minrate = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'MINRATE', & ! tag name
    'MINRATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'minimum shutoff rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_maxrate = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'MAXRATE', & ! tag name
    'MAXRATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'maximum shutoff rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_ratescalingrec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'RATE_SCALINGRECORD', & ! tag name
    'RATESCALINGREC', & ! fortran variable
    'RECORD RATE_SCALING PUMP_ELEVATION SCALING_LENGTH', & ! type
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
    gwfmaw_rate_scaling = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'RATE_SCALING', & ! tag name
    'RATE_SCALING', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'rate scaling', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_pump_elevation = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'PUMP_ELEVATION', & ! tag name
    'PUMP_ELEVATION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'pump elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_scaling_length = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'SCALING_LENGTH', & ! tag name
    'SCALING_LENGTH', & ! fortran variable
    'DOUBLE', & ! type
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
    gwfmaw_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARYRECORD', & ! tag name
    'AUXILIARYRECORD', & ! fortran variable
    'RECORD AUXILIARY AUXNAME AUXVAL', & ! type
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
    gwfmaw_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARY', & ! tag name
    'PERIOD_AUXILIARY', & ! fortran variable
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
    gwfmaw_auxname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'AUXNAME', & ! tag name
    'AUXNAME', & ! fortran variable
    'STRING', & ! type
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
    gwfmaw_auxval = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'AUXVAL', & ! tag name
    'AUXVAL', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'auxiliary variable value', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_maw_param_definitions(*) = &
    [ &
    gwfmaw_auxiliary, &
    gwfmaw_boundnames, &
    gwfmaw_iprpak, &
    gwfmaw_print_head, &
    gwfmaw_iprflow, &
    gwfmaw_ipakcb, &
    gwfmaw_head_filerecord, &
    gwfmaw_head, &
    gwfmaw_headfile, &
    gwfmaw_budfilerec, &
    gwfmaw_budget, &
    gwfmaw_fileout, &
    gwfmaw_budgetfile, &
    gwfmaw_budcsvfilerec, &
    gwfmaw_budgetcsv, &
    gwfmaw_budgetcsvfile, &
    gwfmaw_no_well_storage, &
    gwfmaw_flow_correction, &
    gwfmaw_inonvert, &
    gwfmaw_flowing_wells, &
    gwfmaw_shutdown_theta, &
    gwfmaw_shutdown_kappa, &
    gwfmaw_mfrcsvfilerec, &
    gwfmaw_maw_flw_red_csv, &
    gwfmaw_mfrcsvfile, &
    gwfmaw_ts_filerecord, &
    gwfmaw_ts6, &
    gwfmaw_filein, &
    gwfmaw_ts6_filename, &
    gwfmaw_obs_filerecord, &
    gwfmaw_obs6, &
    gwfmaw_obs6_filename, &
    gwfmaw_mover, &
    gwfmaw_ieffradopt, &
    gwfmaw_nmawwells, &
    gwfmaw_packagedata_ifno, &
    gwfmaw_radius, &
    gwfmaw_bottom, &
    gwfmaw_strt, &
    gwfmaw_condeqn, &
    gwfmaw_ngwfnodes, &
    gwfmaw_aux, &
    gwfmaw_boundname, &
    gwfmaw_conndata_ifno, &
    gwfmaw_icon, &
    gwfmaw_cellid, &
    gwfmaw_scrn_top, &
    gwfmaw_scrn_bot, &
    gwfmaw_hk_skin, &
    gwfmaw_radius_skin, &
    gwfmaw_angledata_ifno, &
    gwfmaw_angledata_icon, &
    gwfmaw_angle, &
    gwfmaw_connlen, &
    gwfmaw_ifno, &
    gwfmaw_status, &
    gwfmaw_flwwellrecord, &
    gwfmaw_flowing_well, &
    gwfmaw_fwelev, &
    gwfmaw_fwcond, &
    gwfmaw_fwrlen, &
    gwfmaw_rate, &
    gwfmaw_well_head, &
    gwfmaw_head_limit, &
    gwfmaw_shutoffrecord, &
    gwfmaw_shut_off, &
    gwfmaw_minrate, &
    gwfmaw_maxrate, &
    gwfmaw_ratescalingrec, &
    gwfmaw_rate_scaling, &
    gwfmaw_pump_elevation, &
    gwfmaw_scaling_length, &
    gwfmaw_auxiliaryrecord, &
    gwfmaw_period_auxiliary, &
    gwfmaw_auxname, &
    gwfmaw_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO RADIUS BOTTOM STRT CONDEQN NGWFNODES AUX '// &
    'BOUNDNAME', & ! type
    'NMAWWELLS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmaw_connectiondata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CONNECTIONDATA', & ! tag name
    'CONNECTIONDATA', & ! fortran variable
    'RECARRAY IFNO ICON CELLID SCRN_TOP SCRN_BOT HK_SKIN '// &
    'RADIUS_SKIN', & ! type
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
    gwfmaw_angledata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'ANGLEDATA', & ! block
    'ANGLEDATA', & ! tag name
    'ANGLEDATA', & ! fortran variable
    'RECARRAY IFNO ICON ANGLE CONN_LENGTH', & ! type
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
    gwfmaw_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY IFNO MAWSETTING', & ! type
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
    gwfmaw_mawsetting = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MAW', & ! subcomponent
    'PERIOD', & ! block
    'MAWSETTING', & ! tag name
    'MAWSETTING', & ! fortran variable
    'KEYSTRING STATUS FLOWING_WELLRECORD RATE WELL_HEAD '// &
    'HEAD_LIMIT SHUTOFFRECORD RATE_SCALINGRECORD AUXILIARYRECORD', & ! type
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
    gwf_maw_aggregate_definitions(*) = &
    [ &
    gwfmaw_packagedata, &
    gwfmaw_connectiondata, &
    gwfmaw_angledata, &
    gwfmaw_perioddata, &
    gwfmaw_mawsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_maw_block_definitions(*) = &
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
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CONNECTIONDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'ANGLEDATA', & ! blockname
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

end module GwfMawInputModule

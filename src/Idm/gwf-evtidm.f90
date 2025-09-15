! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfEvtInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_evt_param_definitions
  public gwf_evt_aggregate_definitions
  public gwf_evt_block_definitions
  public GwfEvtParamFoundType
  public gwf_evt_multi_package
  public gwf_evt_subpackages

  type GwfEvtParamFoundType
    logical :: fixed_cell = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: surfratespec = .false.
    logical :: maxbound = .false.
    logical :: nseg = .false.
    logical :: cellid = .false.
    logical :: surface = .false.
    logical :: rate = .false.
    logical :: depth = .false.
    logical :: pxdp = .false.
    logical :: petm = .false.
    logical :: petm0 = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwfEvtParamFoundType

  logical :: gwf_evt_multi_package = .true.

  character(len=16), parameter :: &
    gwf_evt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfevt_fixed_cell = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'FIXED_CELL', & ! tag name
    'FIXED_CELL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'if cell is dry do not apply evapotranspiration to underlying cell', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'AUXMULTNAME', & ! tag name
    'AUXMULTNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for multiplier', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
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
    gwfevt_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print evapotranspiration rates to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save evapotranspiration rates to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
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
    gwfevt_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
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
    gwfevt_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
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
    gwfevt_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_surfratespec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'OPTIONS', & ! block
    'SURF_RATE_SPECIFIED', & ! tag name
    'SURFRATESPEC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'specify proportion of evapotranspiration rate at ET surface', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of evapotranspiration cells', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_nseg = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'DIMENSIONS', & ! block
    'NSEG', & ! tag name
    'NSEG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of ET segments', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_surface = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'SURFACE', & ! tag name
    'SURFACE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'ET surface', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_rate = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'maximum ET rate', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_depth = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'DEPTH', & ! tag name
    'DEPTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'ET extinction depth', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_pxdp = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'PXDP', & ! tag name
    'PXDP', & ! fortran variable
    'DOUBLE1D', & ! type
    'NSEG-1', & ! shape
    'proportion of ET extinction depth', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_petm = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'PETM', & ! tag name
    'PETM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NSEG-1', & ! shape
    'proportion of maximum ET rate', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_petm0 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'PETM0', & ! tag name
    'PETM0', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'proportion of maximum ET rate at ET surface', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfevt_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'evapotranspiration name', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_evt_param_definitions(*) = &
    [ &
    gwfevt_fixed_cell, &
    gwfevt_auxiliary, &
    gwfevt_auxmultname, &
    gwfevt_boundnames, &
    gwfevt_iprpak, &
    gwfevt_iprflow, &
    gwfevt_ipakcb, &
    gwfevt_ts_filerecord, &
    gwfevt_ts6, &
    gwfevt_filein, &
    gwfevt_ts6_filename, &
    gwfevt_obs_filerecord, &
    gwfevt_obs6, &
    gwfevt_obs6_filename, &
    gwfevt_surfratespec, &
    gwfevt_maxbound, &
    gwfevt_nseg, &
    gwfevt_cellid, &
    gwfevt_surface, &
    gwfevt_rate, &
    gwfevt_depth, &
    gwfevt_pxdp, &
    gwfevt_petm, &
    gwfevt_petm0, &
    gwfevt_auxvar, &
    gwfevt_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfevt_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'EVT', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID SURFACE RATE DEPTH PXDP PETM PETM0 AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_evt_aggregate_definitions(*) = &
    [ &
    gwfevt_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_evt_block_definitions(*) = &
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

end module GwfEvtInputModule

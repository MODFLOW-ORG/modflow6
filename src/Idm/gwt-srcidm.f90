! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtSrcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_src_param_definitions
  public gwt_src_aggregate_definitions
  public gwt_src_block_definitions
  public GwtSrcParamFoundType
  public gwt_src_multi_package
  public gwt_src_subpackages

  type GwtSrcParamFoundType
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
    logical :: highest_sat = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: smassrate = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwtSrcParamFoundType

  logical :: gwt_src_multi_package = .true.

  character(len=16), parameter :: &
    gwt_src_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_auxmultname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'AUXMULTNAME', & ! tag name
    'AUXMULTNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for multiplier', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_boundnames = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_iprpak = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save well flows to budget file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_ts6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'head keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_ts6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_obs6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_obs6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_highest_sat = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'OPTIONS', & ! block
    'HIGHEST_SATURATED', & ! tag name
    'HIGHEST_SAT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'apply source to highest saturated cell', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_maxbound = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of sources', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_cellid = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_smassrate = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'PERIOD', & ! block
    'SMASSRATE', & ! tag name
    'SMASSRATE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'mass source loading rate', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_auxvar = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_boundname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well name', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_src_param_definitions(*) = &
    [ &
    gwtsrc_auxiliary, &
    gwtsrc_auxmultname, &
    gwtsrc_boundnames, &
    gwtsrc_iprpak, &
    gwtsrc_iprflow, &
    gwtsrc_ipakcb, &
    gwtsrc_ts_filerecord, &
    gwtsrc_ts6, &
    gwtsrc_filein, &
    gwtsrc_ts6_filename, &
    gwtsrc_obs_filerecord, &
    gwtsrc_obs6, &
    gwtsrc_obs6_filename, &
    gwtsrc_highest_sat, &
    gwtsrc_maxbound, &
    gwtsrc_cellid, &
    gwtsrc_smassrate, &
    gwtsrc_auxvar, &
    gwtsrc_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtsrc_spd = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SRC', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID SMASSRATE AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_src_aggregate_definitions(*) = &
    [ &
    gwtsrc_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_src_block_definitions(*) = &
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

end module GwtSrcInputModule

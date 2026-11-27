! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfVssbInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_vssb_param_definitions
  public gwf_vssb_aggregate_definitions
  public gwf_vssb_block_definitions
  public GwfVssbParamFoundType
  public gwf_vssb_multi_package
  public gwf_vssb_is_advanced
  public gwf_vssb_subpackages

  type GwfVssbParamFoundType
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: print_flows = .false.
    logical :: ipakcb = .false.
    logical :: pond_filerecord = .false.
    logical :: ponddepth = .false.
    logical :: fileout = .false.
    logical :: ponddepthfile = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: maxbound = .false.
    logical :: packagedata_ifno = .false.
    logical :: cellid = .false.
    logical :: max_pond_depth = .false.
    logical :: boundname = .false.
    logical :: ifno = .false.
    logical :: precipitation = .false.
    logical :: evaporation = .false.
    logical :: aux = .false.
  end type GwfVssbParamFoundType

  logical :: gwf_vssb_multi_package = .true.
  logical :: gwf_vssb_is_advanced = .true.

  character(len=16), parameter :: &
    gwf_vssb_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'AUXMULTNAME', & ! tag name
    'AUXMULTNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for multiplier', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_print_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
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
    gwfvssb_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_pond_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'POND_FILERECORD', & ! tag name
    'POND_FILERECORD', & ! fortran variable
    'RECORD PONDDEPTH FILEOUT PONDDEPTHFILE', & ! type
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
    gwfvssb_ponddepth = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'PONDDEPTH', & ! tag name
    'PONDDEPTH', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'pond depth keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_ponddepthfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'OPTIONS', & ! block
    'PONDDEPTHFILE', & ! tag name
    'PONDDEPTHFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'pond depth file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'DIMENSIONS', & ! block
    'NSURFCELLS', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of surface cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'surface cell id number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
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
    gwfvssb_max_pond_depth = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'MAX_POND_DEPTH', & ! tag name
    'MAX_POND_DEPTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'maximum pond depth', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'surface cell name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'surface cell id number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_precipitation = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PERIOD', & ! block
    'PRECIPITATION', & ! tag name
    'PRECIPITATION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'precipitation rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_evaporation = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PERIOD', & ! block
    'EVAPORATION', & ! tag name
    'EVAPORATION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'potential evaporation rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_aux = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PERIOD', & ! block
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
    gwf_vssb_param_definitions(*) = &
    [ &
    gwfvssb_auxiliary, &
    gwfvssb_auxmultname, &
    gwfvssb_boundnames, &
    gwfvssb_iprpak, &
    gwfvssb_print_flows, &
    gwfvssb_ipakcb, &
    gwfvssb_pond_filerecord, &
    gwfvssb_ponddepth, &
    gwfvssb_fileout, &
    gwfvssb_ponddepthfile, &
    gwfvssb_ts_filerecord, &
    gwfvssb_ts6, &
    gwfvssb_filein, &
    gwfvssb_ts6_filename, &
    gwfvssb_obs_filerecord, &
    gwfvssb_obs6, &
    gwfvssb_obs6_filename, &
    gwfvssb_mover, &
    gwfvssb_maxbound, &
    gwfvssb_packagedata_ifno, &
    gwfvssb_cellid, &
    gwfvssb_max_pond_depth, &
    gwfvssb_boundname, &
    gwfvssb_ifno, &
    gwfvssb_precipitation, &
    gwfvssb_evaporation, &
    gwfvssb_aux &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO CELLID MAX_POND_DEPTH BOUNDNAME', & ! type
    'NSURFCELLS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvssb_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSB', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY IFNO PRECIPITATION EVAPORATION AUX', & ! type
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
    gwf_vssb_aggregate_definitions(*) = &
    [ &
    gwfvssb_packagedata, &
    gwfvssb_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_vssb_block_definitions(*) = &
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
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfVssbInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfCdbInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_cdb_param_definitions
  public chf_cdb_aggregate_definitions
  public chf_cdb_block_definitions
  public ChfCdbParamFoundType
  public chf_cdb_multi_package
  public chf_cdb_subpackages

  type ChfCdbParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: filein = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: idcxs = .false.
    logical :: width = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type ChfCdbParamFoundType

  logical :: chf_cdb_multi_package = .true.

  character(len=16), parameter :: &
    chf_cdb_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfcdb_auxiliary = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_boundnames = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_iprpak = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_iprflow = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_ipakcb = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_filein = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_obs_filerecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_obs6 = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_obs6_filename = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
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
    chfcdb_maxbound = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of critical depth boundaries', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcdb_cellid = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
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
    chfcdb_idcxs = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cross section identifier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcdb_width = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
    'WIDTH', & ! tag name
    'WIDTH', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'width of the zero-depth gradient boundary', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcdb_auxvar = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
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
    chfcdb_boundname = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'zero-depth-gradient boundary name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_cdb_param_definitions(*) = &
    [ &
    chfcdb_auxiliary, &
    chfcdb_boundnames, &
    chfcdb_iprpak, &
    chfcdb_iprflow, &
    chfcdb_ipakcb, &
    chfcdb_filein, &
    chfcdb_obs_filerecord, &
    chfcdb_obs6, &
    chfcdb_obs6_filename, &
    chfcdb_maxbound, &
    chfcdb_cellid, &
    chfcdb_idcxs, &
    chfcdb_width, &
    chfcdb_auxvar, &
    chfcdb_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfcdb_spd = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CDB', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID IDCXS WIDTH AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_cdb_aggregate_definitions(*) = &
    [ &
    chfcdb_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_cdb_block_definitions(*) = &
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

end module ChfCdbInputModule

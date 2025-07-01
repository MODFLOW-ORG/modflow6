! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfGhbgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_ghbg_param_definitions
  public gwf_ghbg_aggregate_definitions
  public gwf_ghbg_block_definitions
  public GwfGhbgParamFoundType
  public gwf_ghbg_multi_package
  public gwf_ghbg_subpackages

  type GwfGhbgParamFoundType
    logical :: readarraygrid = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: filein = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: export_nc = .false.
    logical :: maxbound = .false.
    logical :: bhead = .false.
    logical :: cond = .false.
    logical :: auxvar = .false.
  end type GwfGhbgParamFoundType

  logical :: gwf_ghbg_multi_package = .true.

  character(len=16), parameter :: &
    gwf_ghbg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_readarraygrid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'OPTIONS', & ! block
    'READARRAYGRID', & ! tag name
    'READARRAYGRID', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use array-based grid input', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save GHBG flows to budget file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
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
    gwfghbg_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'OPTIONS', & ! block
    'MOVER', & ! tag name
    'MOVER', & ! fortran variable
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
    gwfghbg_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of general-head boundaries in any stress period', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_bhead = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'PERIOD', & ! block
    'BHEAD', & ! tag name
    'BHEAD', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'boundary head', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_cond = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'PERIOD', & ! block
    'COND', & ! tag name
    'COND', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'boundary conductance', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfghbg_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GHBG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'general-head boundary auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_ghbg_param_definitions(*) = &
    [ &
    gwfghbg_readarraygrid, &
    gwfghbg_auxiliary, &
    gwfghbg_auxmultname, &
    gwfghbg_iprpak, &
    gwfghbg_iprflow, &
    gwfghbg_ipakcb, &
    gwfghbg_obs_filerecord, &
    gwfghbg_obs6, &
    gwfghbg_filein, &
    gwfghbg_obs6_filename, &
    gwfghbg_mover, &
    gwfghbg_export_nc, &
    gwfghbg_maxbound, &
    gwfghbg_bhead, &
    gwfghbg_cond, &
    gwfghbg_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_ghbg_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_ghbg_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
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
    'PERIOD', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfGhbgInputModule

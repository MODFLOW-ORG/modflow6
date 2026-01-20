! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfChdgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_chdg_param_definitions
  public gwf_chdg_aggregate_definitions
  public gwf_chdg_block_definitions
  public GwfChdgParamFoundType
  public gwf_chdg_multi_package
  public gwf_chdg_subpackages

  type GwfChdgParamFoundType
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
    logical :: export_nc = .false.
    logical :: inewton = .false.
    logical :: maxbound = .false.
    logical :: head = .false.
    logical :: auxvar = .false.
  end type GwfChdgParamFoundType

  logical :: gwf_chdg_multi_package = .true.

  character(len=16), parameter :: &
    gwf_chdg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_readarraygrid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'OPTIONS', & ! block
    'READARRAYGRID', & ! tag name
    'READARRAYGRID', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use array-based grid input', & ! longname
    .true., & ! required
    .true., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print CHD flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save CHD flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
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
    gwfchdg_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_inewton = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_NEWTON', & ! tag name
    'INEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'turn off Newton for unconfined cells', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of constant head cells in any stress period', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'PERIOD', & ! block
    'HEAD', & ! tag name
    'HEAD', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'head value assigned to constant head', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfchdg_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CHDG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'constant head auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_chdg_param_definitions(*) = &
    [ &
    gwfchdg_readarraygrid, &
    gwfchdg_auxiliary, &
    gwfchdg_auxmultname, &
    gwfchdg_iprpak, &
    gwfchdg_iprflow, &
    gwfchdg_ipakcb, &
    gwfchdg_obs_filerecord, &
    gwfchdg_obs6, &
    gwfchdg_filein, &
    gwfchdg_obs6_filename, &
    gwfchdg_export_nc, &
    gwfchdg_inewton, &
    gwfchdg_maxbound, &
    gwfchdg_head, &
    gwfchdg_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_chdg_aggregate_definitions(*) = &
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
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_chdg_block_definitions(*) = &
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

end module GwfChdgInputModule

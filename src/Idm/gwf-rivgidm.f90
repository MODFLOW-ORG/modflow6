! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfRivgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_rivg_param_definitions
  public gwf_rivg_aggregate_definitions
  public gwf_rivg_block_definitions
  public GwfRivgParamFoundType
  public gwf_rivg_multi_package
  public gwf_rivg_subpackages

  type GwfRivgParamFoundType
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
    logical :: stage = .false.
    logical :: cond = .false.
    logical :: rbot = .false.
    logical :: auxvar = .false.
  end type GwfRivgParamFoundType

  logical :: gwf_rivg_multi_package = .true.

  character(len=16), parameter :: &
    gwf_rivg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_readarraygrid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save RIV flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
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
    gwfrivg_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of river cells in any stress period', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_stage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'PERIOD', & ! block
    'STAGE', & ! tag name
    'STAGE', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'river stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_cond = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'PERIOD', & ! block
    'COND', & ! tag name
    'COND', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'river conductnace', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_rbot = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'PERIOD', & ! block
    'RBOT', & ! tag name
    'RBOT', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'river bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrivg_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RIVG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'river auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_rivg_param_definitions(*) = &
    [ &
    gwfrivg_readarraygrid, &
    gwfrivg_auxiliary, &
    gwfrivg_auxmultname, &
    gwfrivg_iprpak, &
    gwfrivg_iprflow, &
    gwfrivg_ipakcb, &
    gwfrivg_obs_filerecord, &
    gwfrivg_obs6, &
    gwfrivg_filein, &
    gwfrivg_obs6_filename, &
    gwfrivg_mover, &
    gwfrivg_export_nc, &
    gwfrivg_maxbound, &
    gwfrivg_stage, &
    gwfrivg_cond, &
    gwfrivg_rbot, &
    gwfrivg_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_rivg_aggregate_definitions(*) = &
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
    gwf_rivg_block_definitions(*) = &
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

end module GwfRivgInputModule

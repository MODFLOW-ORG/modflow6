! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfRchaInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_rcha_param_definitions
  public gwf_rcha_aggregate_definitions
  public gwf_rcha_block_definitions
  public GwfRchaParamFoundType
  public gwf_rcha_multi_package
  public gwf_rcha_subpackages

  type GwfRchaParamFoundType
    logical :: readasarrays = .false.
    logical :: fixed_cell = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: tas_filerecord = .false.
    logical :: tas6 = .false.
    logical :: filein = .false.
    logical :: tas6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: export_nc = .false.
    logical :: irch = .false.
    logical :: recharge = .false.
    logical :: auxvar = .false.
  end type GwfRchaParamFoundType

  logical :: gwf_rcha_multi_package = .true.

  character(len=16), parameter :: &
    gwf_rcha_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_readasarrays = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'READASARRAYS', & ! tag name
    'READASARRAYS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use array-based input', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_fixed_cell = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'FIXED_CELL', & ! tag name
    'FIXED_CELL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'if cell is dry do not apply recharge to underlying cell', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print recharge rates to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save CHD flows to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_tas_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS_FILERECORD', & ! tag name
    'TAS_FILERECORD', & ! fortran variable
    'RECORD TAS6 FILEIN TAS6_FILENAME', & ! type
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
    gwfrcha_tas6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6', & ! tag name
    'TAS6', & ! fortran variable
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
    gwfrcha_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_tas6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6_FILENAME', & ! tag name
    'TAS6_FILENAME', & ! fortran variable
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
    gwfrcha_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
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
    gwfrcha_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_irch = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'PERIOD', & ! block
    'IRCH', & ! tag name
    'IRCH', & ! fortran variable
    'INTEGER1D', & ! type
    'NCPL', & ! shape
    'layer number for recharge', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_recharge = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'PERIOD', & ! block
    'RECHARGE', & ! tag name
    'RECHARGE', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCPL', & ! shape
    'recharge rate', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfrcha_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'RCHA', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NCPL', & ! shape
    'recharge auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_rcha_param_definitions(*) = &
    [ &
    gwfrcha_readasarrays, &
    gwfrcha_fixed_cell, &
    gwfrcha_auxiliary, &
    gwfrcha_auxmultname, &
    gwfrcha_iprpak, &
    gwfrcha_iprflow, &
    gwfrcha_ipakcb, &
    gwfrcha_tas_filerecord, &
    gwfrcha_tas6, &
    gwfrcha_filein, &
    gwfrcha_tas6_filename, &
    gwfrcha_obs_filerecord, &
    gwfrcha_obs6, &
    gwfrcha_obs6_filename, &
    gwfrcha_export_nc, &
    gwfrcha_irch, &
    gwfrcha_recharge, &
    gwfrcha_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_rcha_aggregate_definitions(*) = &
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
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_rcha_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
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

end module GwfRchaInputModule

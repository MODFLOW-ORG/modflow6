! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfDrngInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_drng_param_definitions
  public gwf_drng_aggregate_definitions
  public gwf_drng_block_definitions
  public GwfDrngParamFoundType
  public gwf_drng_multi_package
  public gwf_drng_subpackages

  type GwfDrngParamFoundType
    logical :: readarraygrid = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: auxdepthname = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: filein = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: export_nc = .false.
    logical :: icubicsfac = .false.
    logical :: maxbound = .false.
    logical :: elev = .false.
    logical :: cond = .false.
    logical :: auxvar = .false.
  end type GwfDrngParamFoundType

  logical :: gwf_drng_multi_package = .true.

  character(len=16), parameter :: &
    gwf_drng_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_readarraygrid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_auxdepthname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'OPTIONS', & ! block
    'AUXDEPTHNAME', & ! tag name
    'AUXDEPTHNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for drainage depth', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save DRNG flows to budget file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
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
    gwfdrng_icubicsfac = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_CUBIC_SCALING', & ! tag name
    'ICUBICSFAC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'cubic-scaling', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of drain cells in any stress period', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_elev = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'PERIOD', & ! block
    'ELEV', & ! tag name
    'ELEV', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'drain elevation', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_cond = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'PERIOD', & ! block
    'COND', & ! tag name
    'COND', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'drain conductance', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdrng_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DRNG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'drain auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_drng_param_definitions(*) = &
    [ &
    gwfdrng_readarraygrid, &
    gwfdrng_auxiliary, &
    gwfdrng_auxmultname, &
    gwfdrng_auxdepthname, &
    gwfdrng_iprpak, &
    gwfdrng_iprflow, &
    gwfdrng_ipakcb, &
    gwfdrng_obs_filerecord, &
    gwfdrng_obs6, &
    gwfdrng_filein, &
    gwfdrng_obs6_filename, &
    gwfdrng_mover, &
    gwfdrng_export_nc, &
    gwfdrng_icubicsfac, &
    gwfdrng_maxbound, &
    gwfdrng_elev, &
    gwfdrng_cond, &
    gwfdrng_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_drng_aggregate_definitions(*) = &
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
    gwf_drng_block_definitions(*) = &
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

end module GwfDrngInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfWelgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_welg_param_definitions
  public gwf_welg_aggregate_definitions
  public gwf_welg_block_definitions
  public GwfWelgParamFoundType
  public gwf_welg_multi_package
  public gwf_welg_subpackages

  type GwfWelgParamFoundType
    logical :: readarraygrid = .false.
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: flowred = .false.
    logical :: afrcsv_rec = .false.
    logical :: afrcsv = .false.
    logical :: fileout = .false.
    logical :: afrcsvfile = .false.
    logical :: obs_filerecord = .false.
    logical :: filein = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: export_nc = .false.
    logical :: maxbound = .false.
    logical :: q = .false.
    logical :: auxvar = .false.
  end type GwfWelgParamFoundType

  logical :: gwf_welg_multi_package = .true.

  character(len=16), parameter :: &
    gwf_welg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_readarraygrid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_flowred = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'OPTIONS', & ! block
    'AUTO_FLOW_REDUCE', & ! tag name
    'FLOWRED', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'cell fractional thickness for reduced pumping', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_afrcsv_rec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'OPTIONS', & ! block
    'AFRCSV_FILERECORD', & ! tag name
    'AFRCSV_REC', & ! fortran variable
    'RECORD AUTO_FLOW_REDUCE_CSV FILEOUT AFRCSVFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_afrcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'OPTIONS', & ! block
    'AUTO_FLOW_REDUCE_CSV', & ! tag name
    'AFRCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
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
    gwfwelg_afrcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'OPTIONS', & ! block
    'AFRCSVFILE', & ! tag name
    'AFRCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
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
    gwfwelg_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of wells in any stress period', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_q = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'PERIOD', & ! block
    'Q', & ! tag name
    'Q', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'well rate', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfwelg_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'WELG', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'well auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_welg_param_definitions(*) = &
    [ &
    gwfwelg_readarraygrid, &
    gwfwelg_auxiliary, &
    gwfwelg_auxmultname, &
    gwfwelg_iprpak, &
    gwfwelg_iprflow, &
    gwfwelg_ipakcb, &
    gwfwelg_flowred, &
    gwfwelg_afrcsv_rec, &
    gwfwelg_afrcsv, &
    gwfwelg_fileout, &
    gwfwelg_afrcsvfile, &
    gwfwelg_obs_filerecord, &
    gwfwelg_filein, &
    gwfwelg_obs6, &
    gwfwelg_obs6_filename, &
    gwfwelg_mover, &
    gwfwelg_export_nc, &
    gwfwelg_maxbound, &
    gwfwelg_q, &
    gwfwelg_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_welg_aggregate_definitions(*) = &
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
    gwf_welg_block_definitions(*) = &
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

end module GwfWelgInputModule

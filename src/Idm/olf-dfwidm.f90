! ** Do Not Modify! MODFLOW 6 system generated file. **
module OlfDfwInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public olf_dfw_param_definitions
  public olf_dfw_aggregate_definitions
  public olf_dfw_block_definitions
  public OlfDfwParamFoundType
  public olf_dfw_multi_package
  public olf_dfw_subpackages

  type OlfDfwParamFoundType
    logical :: icentral = .false.
    logical :: lengthconv = .false.
    logical :: timeconv = .false.
    logical :: ipakcb = .false.
    logical :: iprflow = .false.
    logical :: isavvelocity = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: filein = .false.
    logical :: obs6_filename = .false.
    logical :: export_ascii = .false.
    logical :: iswrcond = .false.
    logical :: manningsn = .false.
    logical :: idcxs = .false.
  end type OlfDfwParamFoundType

  logical :: olf_dfw_multi_package = .false.

  character(len=16), parameter :: &
    olf_dfw_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    olfdfw_icentral = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'CENTRAL_IN_SPACE', & ! tag name
    'ICENTRAL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use central in space weighting', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_lengthconv = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_CONVERSION', & ! tag name
    'LENGTHCONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'length conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_timeconv = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_CONVERSION', & ! tag name
    'TIMECONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'time conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_ipakcb = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save DFW flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_iprflow = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print DFW flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_isavvelocity = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_VELOCITY', & ! tag name
    'ISAVVELOCITY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save velocity', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_obs_filerecord = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
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
    olfdfw_obs6 = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
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
    olfdfw_filein = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
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
    olfdfw_obs6_filename = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
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
    olfdfw_export_ascii = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_iswrcond = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_SWR_CONDUCTANCE', & ! tag name
    'ISWRCOND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use SWR conductance formulation', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_manningsn = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'MANNINGSN', & ! tag name
    'MANNINGSN', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfdfw_idcxs = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'cross section number', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olf_dfw_param_definitions(*) = &
    [ &
    olfdfw_icentral, &
    olfdfw_lengthconv, &
    olfdfw_timeconv, &
    olfdfw_ipakcb, &
    olfdfw_iprflow, &
    olfdfw_isavvelocity, &
    olfdfw_obs_filerecord, &
    olfdfw_obs6, &
    olfdfw_filein, &
    olfdfw_obs6_filename, &
    olfdfw_export_ascii, &
    olfdfw_iswrcond, &
    olfdfw_manningsn, &
    olfdfw_idcxs &
    ]

  type(InputParamDefinitionType), parameter :: &
    olf_dfw_aggregate_definitions(*) = &
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
    olf_dfw_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module OlfDfwInputModule

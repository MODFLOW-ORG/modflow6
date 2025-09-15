! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfDfwInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_dfw_param_definitions
  public chf_dfw_aggregate_definitions
  public chf_dfw_block_definitions
  public ChfDfwParamFoundType
  public chf_dfw_multi_package
  public chf_dfw_subpackages

  type ChfDfwParamFoundType
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
  end type ChfDfwParamFoundType

  logical :: chf_dfw_multi_package = .false.

  character(len=16), parameter :: &
    chf_dfw_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfdfw_icentral = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'CENTRAL_IN_SPACE', & ! tag name
    'ICENTRAL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use central in space weighting', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_lengthconv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_CONVERSION', & ! tag name
    'LENGTHCONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'length conversion factor', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_timeconv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_CONVERSION', & ! tag name
    'TIMECONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'time conversion factor', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_ipakcb = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save DFW flows', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_iprflow = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print DFW flows to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_isavvelocity = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_VELOCITY', & ! tag name
    'ISAVVELOCITY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save velocity', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_obs_filerecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
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
    chfdfw_obs6 = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
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
    chfdfw_filein = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
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
    chfdfw_obs6_filename = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
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
    chfdfw_export_ascii = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_iswrcond = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_SWR_CONDUCTANCE', & ! tag name
    'ISWRCOND', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use SWR conductance formulation', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_manningsn = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'MANNINGSN', & ! tag name
    'MANNINGSN', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdfw_idcxs = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DFW', & ! subcomponent
    'GRIDDATA', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'cross section number', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_dfw_param_definitions(*) = &
    [ &
    chfdfw_icentral, &
    chfdfw_lengthconv, &
    chfdfw_timeconv, &
    chfdfw_ipakcb, &
    chfdfw_iprflow, &
    chfdfw_isavvelocity, &
    chfdfw_obs_filerecord, &
    chfdfw_obs6, &
    chfdfw_filein, &
    chfdfw_obs6_filename, &
    chfdfw_export_ascii, &
    chfdfw_iswrcond, &
    chfdfw_manningsn, &
    chfdfw_idcxs &
    ]

  type(InputParamDefinitionType), parameter :: &
    chf_dfw_aggregate_definitions(*) = &
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
    chf_dfw_block_definitions(*) = &
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

end module ChfDfwInputModule

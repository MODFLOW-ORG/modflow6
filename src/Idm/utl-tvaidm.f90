! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlTvaInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_tva_param_definitions
  public utl_tva_aggregate_definitions
  public utl_tva_block_definitions
  public UtlTvaParamFoundType
  public utl_tva_multi_package
  public utl_tva_is_advanced
  public utl_tva_subpackages

  type UtlTvaParamFoundType
    logical :: readasarrays = .false.
    logical :: auxiliary = .false.
    logical :: iprpak = .false.
    logical :: tas_filerecord = .false.
    logical :: tas6 = .false.
    logical :: filein = .false.
    logical :: tas6_filename = .false.
    logical :: export_nc = .false.
    logical :: auxvar = .false.
  end type UtlTvaParamFoundType

  logical :: utl_tva_multi_package = .true.
  logical :: utl_tva_is_advanced = .false.

  character(len=16), parameter :: &
    utl_tva_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utltva_readasarrays = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
    'OPTIONS', & ! block
    'READASARRAYS', & ! tag name
    'READASARRAYS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use array-based input', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utltva_auxiliary = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
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
    utltva_iprpak = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
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
    utltva_tas_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS_FILERECORD', & ! tag name
    'TAS_FILERECORD', & ! fortran variable
    'RECORD TAS6 FILEIN TAS6_FILENAME', & ! type
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
    utltva_tas6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6', & ! tag name
    'TAS6', & ! fortran variable
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
    utltva_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
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
    utltva_tas6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
    'OPTIONS', & ! block
    'TAS6_FILENAME', & ! tag name
    'TAS6_FILENAME', & ! fortran variable
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
    utltva_export_nc = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
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
    utltva_auxvar = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'TVA', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE2D', & ! type
    'NAUX NODES', & ! shape
    'auxiliary variable iaux', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_tva_param_definitions(*) = &
    [ &
    utltva_readasarrays, &
    utltva_auxiliary, &
    utltva_iprpak, &
    utltva_tas_filerecord, &
    utltva_tas6, &
    utltva_filein, &
    utltva_tas6_filename, &
    utltva_export_nc, &
    utltva_auxvar &
    ]

  type(InputParamDefinitionType), parameter :: &
    utl_tva_aggregate_definitions(*) = &
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
    utl_tva_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module UtlTvaInputModule

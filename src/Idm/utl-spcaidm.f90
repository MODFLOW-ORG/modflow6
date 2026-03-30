! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlSpcaInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_spca_param_definitions
  public utl_spca_aggregate_definitions
  public utl_spca_block_definitions
  public UtlSpcaParamFoundType
  public utl_spca_multi_package
  public utl_spca_subpackages

  type UtlSpcaParamFoundType
    logical :: readasarrays = .false.
    logical :: print_input = .false.
    logical :: tas_filerecord = .false.
    logical :: tas6 = .false.
    logical :: filein = .false.
    logical :: tas6_filename = .false.
    logical :: concentration = .false.
    logical :: temperature = .false.
  end type UtlSpcaParamFoundType

  logical :: utl_spca_multi_package = .true.

  character(len=16), parameter :: &
    utl_spca_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlspca_readasarrays = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
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
    utlspca_print_input = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
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
    utlspca_tas_filerecord = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
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
    utlspca_tas6 = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
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
    utlspca_filein = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
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
    utlspca_tas6_filename = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
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
    utlspca_concentration = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCPL', & ! shape
    'concentration', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlspca_temperature = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SPCA', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCPL', & ! shape
    'temperature', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_spca_param_definitions(*) = &
    [ &
    utlspca_readasarrays, &
    utlspca_print_input, &
    utlspca_tas_filerecord, &
    utlspca_tas6, &
    utlspca_filein, &
    utlspca_tas6_filename, &
    utlspca_concentration, &
    utlspca_temperature &
    ]

  type(InputParamDefinitionType), parameter :: &
    utl_spca_aggregate_definitions(*) = &
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
    utl_spca_block_definitions(*) = &
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

end module UtlSpcaInputModule

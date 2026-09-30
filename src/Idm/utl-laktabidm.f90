! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlLaktabInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_laktab_param_definitions
  public utl_laktab_aggregate_definitions
  public utl_laktab_block_definitions
  public UtlLaktabParamFoundType
  public utl_laktab_multi_package
  public utl_laktab_is_advanced
  public utl_laktab_subpackages

  type UtlLaktabParamFoundType
    logical :: nrow = .false.
    logical :: ncol = .false.
    logical :: stage = .false.
    logical :: volume = .false.
    logical :: sarea = .false.
    logical :: barea = .false.
  end type UtlLaktabParamFoundType

  logical :: utl_laktab_multi_package = .true.
  logical :: utl_laktab_is_advanced = .false.

  character(len=16), parameter :: &
    utl_laktab_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utllaktab_nrow = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'DIMENSIONS', & ! block
    'NROW', & ! tag name
    'NROW', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of table rows', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utllaktab_ncol = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'DIMENSIONS', & ! block
    'NCOL', & ! tag name
    'NCOL', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of table columns', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utllaktab_stage = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'TABLE', & ! block
    'STAGE', & ! tag name
    'STAGE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'lake stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utllaktab_volume = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'TABLE', & ! block
    'VOLUME', & ! tag name
    'VOLUME', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'lake volume', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utllaktab_sarea = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'TABLE', & ! block
    'SAREA', & ! tag name
    'SAREA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'lake surface area', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utllaktab_barea = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'TABLE', & ! block
    'BAREA', & ! tag name
    'BAREA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'lake-GWF exchange area', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_laktab_param_definitions(*) = &
    [ &
    utllaktab_nrow, &
    utllaktab_ncol, &
    utllaktab_stage, &
    utllaktab_volume, &
    utllaktab_sarea, &
    utllaktab_barea &
    ]

  type(InputParamDefinitionType), parameter :: &
    utllaktab_table = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'LAKTAB', & ! subcomponent
    'TABLE', & ! block
    'TABLE', & ! tag name
    'TABLE', & ! fortran variable
    'RECARRAY STAGE VOLUME SAREA BAREA', & ! type
    'NROW', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_laktab_aggregate_definitions(*) = &
    [ &
    utllaktab_table &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_laktab_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'TABLE', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module UtlLaktabInputModule

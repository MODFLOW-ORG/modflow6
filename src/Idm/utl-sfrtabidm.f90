! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlSfrtabInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_sfrtab_param_definitions
  public utl_sfrtab_aggregate_definitions
  public utl_sfrtab_block_definitions
  public UtlSfrtabParamFoundType
  public utl_sfrtab_multi_package
  public utl_sfrtab_is_advanced
  public utl_sfrtab_subpackages

  type UtlSfrtabParamFoundType
    logical :: nrow = .false.
    logical :: ncol = .false.
    logical :: xfraction = .false.
    logical :: height = .false.
    logical :: manfraction = .false.
  end type UtlSfrtabParamFoundType

  logical :: utl_sfrtab_multi_package = .true.
  logical :: utl_sfrtab_is_advanced = .false.

  character(len=16), parameter :: &
    utl_sfrtab_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlsfrtab_nrow = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
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
    utlsfrtab_ncol = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
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
    utlsfrtab_xfraction = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
    'TABLE', & ! block
    'XFRACTION', & ! tag name
    'XFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fractional width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlsfrtab_height = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
    'TABLE', & ! block
    'HEIGHT', & ! tag name
    'HEIGHT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'depth', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlsfrtab_manfraction = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
    'TABLE', & ! block
    'MANFRACTION', & ! tag name
    'MANFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'Mannings roughness coefficient', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_sfrtab_param_definitions(*) = &
    [ &
    utlsfrtab_nrow, &
    utlsfrtab_ncol, &
    utlsfrtab_xfraction, &
    utlsfrtab_height, &
    utlsfrtab_manfraction &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlsfrtab_table = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'SFRTAB', & ! subcomponent
    'TABLE', & ! block
    'TABLE', & ! tag name
    'TABLE', & ! fortran variable
    'RECARRAY XFRACTION HEIGHT MANFRACTION', & ! type
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
    utl_sfrtab_aggregate_definitions(*) = &
    [ &
    utlsfrtab_table &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_sfrtab_block_definitions(*) = &
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

end module UtlSfrtabInputModule

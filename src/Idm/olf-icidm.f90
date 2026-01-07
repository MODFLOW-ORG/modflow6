! ** Do Not Modify! MODFLOW 6 system generated file. **
module OlfIcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public olf_ic_param_definitions
  public olf_ic_aggregate_definitions
  public olf_ic_block_definitions
  public OlfIcParamFoundType
  public olf_ic_multi_package
  public olf_ic_subpackages

  type OlfIcParamFoundType
    logical :: export_ascii = .false.
    logical :: strt = .false.
  end type OlfIcParamFoundType

  logical :: olf_ic_multi_package = .false.

  character(len=16), parameter :: &
    olf_ic_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    olfic_export_ascii = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'IC', & ! subcomponent
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
    olfic_strt = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'IC', & ! subcomponent
    'GRIDDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'starting concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olf_ic_param_definitions(*) = &
    [ &
    olfic_export_ascii, &
    olfic_strt &
    ]

  type(InputParamDefinitionType), parameter :: &
    olf_ic_aggregate_definitions(*) = &
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
    olf_ic_block_definitions(*) = &
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

end module OlfIcInputModule

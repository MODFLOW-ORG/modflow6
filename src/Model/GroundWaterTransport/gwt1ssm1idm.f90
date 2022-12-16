module GwtSsmInputModule
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_ssm_param_definitions
  public gwt_ssm_aggregate_definitions
  public gwt_ssm_block_definitions
  public GwtSsmParamFoundType

  type GwtSsmParamFoundType
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: pname_sources = .false.
    logical :: srctype = .false.
    logical :: auxname = .false.
    logical :: pname_fileinput = .false.
    logical :: spc6 = .false.
    logical :: filein = .false.
    logical :: spc6_filename = .false.
    logical :: mixed = .false.
  end type GwtSsmParamFoundType

  type(InputParamDefinitionType), parameter :: &
    gwtssm_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_pname_sources = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'PNAME', & ! tag name
    'PNAME_SOURCES', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_srctype = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'SRCTYPE', & ! tag name
    'SRCTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_auxname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'AUXNAME', & ! tag name
    'AUXNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_pname_fileinput = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'PNAME', & ! tag name
    'PNAME_FILEINPUT', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_spc6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'SPC6', & ! tag name
    'SPC6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_spc6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'SPC6_FILENAME', & ! tag name
    'SPC6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_mixed = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'MIXED', & ! tag name
    'MIXED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_ssm_param_definitions(*) = &
    [ &
    gwtssm_iprflow, &
    gwtssm_ipakcb, &
    gwtssm_pname_sources, &
    gwtssm_srctype, &
    gwtssm_auxname, &
    gwtssm_pname_fileinput, &
    gwtssm_spc6, &
    gwtssm_filein, &
    gwtssm_spc6_filename, &
    gwtssm_mixed &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtssm_sources = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'SOURCES', & ! tag name
    'SOURCES', & ! fortran variable
    'RECARRAY PNAME SRCTYPE AUXNAME', & ! type
    '', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_fileinput = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'FILEINPUT', & ! tag name
    'FILEINPUT', & ! fortran variable
    'RECARRAY PNAME SPC6 FILEIN SPC6_FILENAME MIXED', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false. & ! layered
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_ssm_aggregate_definitions(*) = &
    [ &
    gwtssm_sources, &
    gwtssm_fileinput &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_ssm_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false. & ! aggregate
    ), &
    InputBlockDefinitionType( &
    'SOURCES', & ! blockname
    .true., & ! required
    .true. & ! aggregate
    ), &
    InputBlockDefinitionType( &
    'FILEINPUT', & ! blockname
    .false., & ! required
    .true. & ! aggregate
    ) &
    ]

end module GwtSsmInputModule

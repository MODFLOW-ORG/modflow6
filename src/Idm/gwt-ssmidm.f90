! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtSsmInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_ssm_param_definitions
  public gwt_ssm_aggregate_definitions
  public gwt_ssm_block_definitions
  public GwtSsmParamFoundType
  public gwt_ssm_multi_package
  public gwt_ssm_subpackages

  type GwtSsmParamFoundType
    logical :: print_flows = .false.
    logical :: save_flows = .false.
    logical :: pname_sources = .false.
    logical :: srctype = .false.
    logical :: auxname = .false.
    logical :: pname = .false.
    logical :: spc6 = .false.
    logical :: filein = .false.
    logical :: spc6_filename = .false.
    logical :: mixed = .false.
  end type GwtSsmParamFoundType

  logical :: gwt_ssm_multi_package = .false.

  character(len=16), parameter :: &
    gwt_ssm_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtssm_print_flows = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_save_flows = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flows to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'package name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'source type', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'auxiliary variable name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtssm_pname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'head keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'spc file name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    'mixed keyword', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_ssm_param_definitions(*) = &
    [ &
    gwtssm_print_flows, &
    gwtssm_save_flows, &
    gwtssm_pname_sources, &
    gwtssm_srctype, &
    gwtssm_auxname, &
    gwtssm_pname, &
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
    'package list', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
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
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'SOURCES', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'FILEINPUT', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwtSsmInputModule

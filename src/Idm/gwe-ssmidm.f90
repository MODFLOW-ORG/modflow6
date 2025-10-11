! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweSsmInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_ssm_param_definitions
  public gwe_ssm_aggregate_definitions
  public gwe_ssm_block_definitions
  public GweSsmParamFoundType
  public gwe_ssm_multi_package
  public gwe_ssm_subpackages

  type GweSsmParamFoundType
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
  end type GweSsmParamFoundType

  logical :: gwe_ssm_multi_package = .false.

  character(len=16), parameter :: &
    gwe_ssm_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwessm_print_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_save_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_pname_sources = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'PNAME', & ! tag name
    'PNAME_SOURCES', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_srctype = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'SRCTYPE', & ! tag name
    'SRCTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'source type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_auxname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'AUXNAME', & ! tag name
    'AUXNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'auxiliary variable name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_pname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_spc6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'SPC6', & ! tag name
    'SPC6', & ! fortran variable
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
    gwessm_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
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
    gwessm_spc6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'SPC6_FILENAME', & ! tag name
    'SPC6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'spc file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_mixed = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'MIXED', & ! tag name
    'MIXED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'mixed keyword', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_ssm_param_definitions(*) = &
    [ &
    gwessm_print_flows, &
    gwessm_save_flows, &
    gwessm_pname_sources, &
    gwessm_srctype, &
    gwessm_auxname, &
    gwessm_pname, &
    gwessm_spc6, &
    gwessm_filein, &
    gwessm_spc6_filename, &
    gwessm_mixed &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwessm_sources = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'SOURCES', & ! block
    'SOURCES', & ! tag name
    'SOURCES', & ! fortran variable
    'RECARRAY PNAME SRCTYPE AUXNAME', & ! type
    '', & ! shape
    'package list', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwessm_fileinput = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SSM', & ! subcomponent
    'FILEINPUT', & ! block
    'FILEINPUT', & ! tag name
    'FILEINPUT', & ! fortran variable
    'RECARRAY PNAME SPC6 FILEIN SPC6_FILENAME MIXED', & ! type
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
    gwe_ssm_aggregate_definitions(*) = &
    [ &
    gwessm_sources, &
    gwessm_fileinput &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_ssm_block_definitions(*) = &
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

end module GweSsmInputModule

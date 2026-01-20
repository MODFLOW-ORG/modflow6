! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfNamInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_nam_param_definitions
  public chf_nam_aggregate_definitions
  public chf_nam_block_definitions
  public ChfNamParamFoundType
  public chf_nam_multi_package
  public chf_nam_subpackages

  type ChfNamParamFoundType
    logical :: list = .false.
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: save_flows = .false.
    logical :: newtonoptions = .false.
    logical :: newton = .false.
    logical :: under_relaxation = .false.
    logical :: ftype = .false.
    logical :: fname = .false.
    logical :: pname = .false.
  end type ChfNamParamFoundType

  logical :: chf_nam_multi_package = .false.

  character(len=16), parameter :: &
    chf_nam_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfnam_list = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'LIST', & ! tag name
    'LIST', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_print_input = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
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
    chfnam_print_flows = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
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
    chfnam_save_flows = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save flows for all packages to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_newtonoptions = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTONOPTIONS', & ! tag name
    'NEWTONOPTIONS', & ! fortran variable
    'RECORD NEWTON UNDER_RELAXATION', & ! type
    '', & ! shape
    'newton keyword and options', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_newton = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTON', & ! tag name
    'NEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson formulation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_under_relaxation = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'UNDER_RELAXATION', & ! tag name
    'UNDER_RELAXATION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson UNDER_RELAXATION option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_ftype = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FTYPE', & ! tag name
    'FTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_fname = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FNAME', & ! tag name
    'FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfnam_pname = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'user name for package', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_nam_param_definitions(*) = &
    [ &
    chfnam_list, &
    chfnam_print_input, &
    chfnam_print_flows, &
    chfnam_save_flows, &
    chfnam_newtonoptions, &
    chfnam_newton, &
    chfnam_under_relaxation, &
    chfnam_ftype, &
    chfnam_fname, &
    chfnam_pname &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfnam_packages = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PACKAGES', & ! tag name
    'PACKAGES', & ! fortran variable
    'RECARRAY FTYPE FNAME PNAME', & ! type
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
    chf_nam_aggregate_definitions(*) = &
    [ &
    chfnam_packages &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_nam_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGES', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ChfNamInputModule

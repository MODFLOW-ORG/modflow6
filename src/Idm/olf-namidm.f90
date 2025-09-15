! ** Do Not Modify! MODFLOW 6 system generated file. **
module OlfNamInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public olf_nam_param_definitions
  public olf_nam_aggregate_definitions
  public olf_nam_block_definitions
  public OlfNamParamFoundType
  public olf_nam_multi_package
  public olf_nam_subpackages

  type OlfNamParamFoundType
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
  end type OlfNamParamFoundType

  logical :: olf_nam_multi_package = .false.

  character(len=16), parameter :: &
    olf_nam_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    olfnam_list = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'LIST', & ! tag name
    'LIST', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_print_input = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_print_flows = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
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
    olfnam_save_flows = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save flows for all packages to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_newtonoptions = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTONOPTIONS', & ! tag name
    'NEWTONOPTIONS', & ! fortran variable
    'RECORD NEWTON UNDER_RELAXATION', & ! type
    '', & ! shape
    'newton keyword and options', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_newton = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTON', & ! tag name
    'NEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson formulation', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_under_relaxation = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'UNDER_RELAXATION', & ! tag name
    'UNDER_RELAXATION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson UNDER_RELAXATION option', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_ftype = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FTYPE', & ! tag name
    'FTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package type', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_fname = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FNAME', & ! tag name
    'FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olfnam_pname = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'user name for package', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    olf_nam_param_definitions(*) = &
    [ &
    olfnam_list, &
    olfnam_print_input, &
    olfnam_print_flows, &
    olfnam_save_flows, &
    olfnam_newtonoptions, &
    olfnam_newton, &
    olfnam_under_relaxation, &
    olfnam_ftype, &
    olfnam_fname, &
    olfnam_pname &
    ]

  type(InputParamDefinitionType), parameter :: &
    olfnam_packages = InputParamDefinitionType &
    ( &
    'OLF', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PACKAGES', & ! tag name
    'PACKAGES', & ! fortran variable
    'RECARRAY FTYPE FNAME PNAME', & ! type
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
    olf_nam_aggregate_definitions(*) = &
    [ &
    olfnam_packages &
    ]

  type(InputBlockDefinitionType), parameter :: &
    olf_nam_block_definitions(*) = &
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

end module OlfNamInputModule

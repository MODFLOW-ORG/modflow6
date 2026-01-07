! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfOcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_oc_param_definitions
  public gwf_oc_aggregate_definitions
  public gwf_oc_block_definitions
  public GwfOcParamFoundType
  public gwf_oc_multi_package
  public gwf_oc_subpackages

  type GwfOcParamFoundType
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: headfilerec = .false.
    logical :: head = .false.
    logical :: headfile = .false.
    logical :: headprintrec = .false.
    logical :: print_format = .false.
    logical :: formatrecord = .false.
    logical :: columns = .false.
    logical :: width = .false.
    logical :: digits = .false.
    logical :: format = .false.
    logical :: saverecord = .false.
    logical :: save = .false.
    logical :: printrecord = .false.
    logical :: print = .false.
    logical :: rtype = .false.
    logical :: ocsetting = .false.
    logical :: all = .false.
    logical :: first = .false.
    logical :: last = .false.
    logical :: frequency = .false.
    logical :: steps = .false.
  end type GwfOcParamFoundType

  logical :: gwf_oc_multi_package = .false.

  character(len=16), parameter :: &
    gwf_oc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfoc_budfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET_FILERECORD', & ! tag name
    'BUDFILEREC', & ! fortran variable
    'RECORD BUDGET FILEOUT BUDGETFILE', & ! type
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
    gwfoc_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET', & ! tag name
    'BUDGET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
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
    gwfoc_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETFILE', & ! tag name
    'BUDGETFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV_FILERECORD', & ! tag name
    'BUDCSVFILEREC', & ! fortran variable
    'RECORD BUDGETCSV FILEOUT BUDGETCSVFILE', & ! type
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
    gwfoc_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV', & ! tag name
    'BUDGETCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSVFILE', & ! tag name
    'BUDGETCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_headfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'HEAD_FILERECORD', & ! tag name
    'HEADFILEREC', & ! fortran variable
    'RECORD HEAD FILEOUT HEADFILE', & ! type
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
    gwfoc_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'HEAD', & ! tag name
    'HEAD', & ! fortran variable
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
    gwfoc_headfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'HEADFILE', & ! tag name
    'HEADFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_headprintrec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'HEADPRINTRECORD', & ! tag name
    'HEADPRINTREC', & ! fortran variable
    'RECORD HEAD PRINT_FORMAT FORMATRECORD', & ! type
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
    gwfoc_print_format = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FORMAT', & ! tag name
    'PRINT_FORMAT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that a print format follows', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_formatrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'FORMATRECORD', & ! tag name
    'FORMATRECORD', & ! fortran variable
    'RECORD COLUMNS WIDTH DIGITS FORMAT', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_columns = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'COLUMNS', & ! tag name
    'COLUMNS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_width = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'WIDTH', & ! tag name
    'WIDTH', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'width for each number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_digits = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'DIGITS', & ! tag name
    'DIGITS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of digits', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_format = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'FORMAT', & ! tag name
    'FORMAT', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'write format', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_saverecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'SAVERECORD', & ! tag name
    'SAVERECORD', & ! fortran variable
    'RECORD SAVE RTYPE OCSETTING', & ! type
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
    gwfoc_save = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'SAVE', & ! tag name
    'SAVE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_printrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'PRINTRECORD', & ! tag name
    'PRINTRECORD', & ! fortran variable
    'RECORD PRINT RTYPE OCSETTING', & ! type
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
    gwfoc_print = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'PRINT', & ! tag name
    'PRINT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_rtype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'RTYPE', & ! tag name
    'RTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'record type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_ocsetting = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'OCSETTING', & ! tag name
    'OCSETTING', & ! fortran variable
    'KEYSTRING ALL FIRST LAST FREQUENCY STEPS', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_all = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'ALL', & ! tag name
    'ALL', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_first = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'FIRST', & ! tag name
    'FIRST', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_last = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'LAST', & ! tag name
    'LAST', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_frequency = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'FREQUENCY', & ! tag name
    'FREQUENCY', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfoc_steps = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'OC', & ! subcomponent
    'PERIOD', & ! block
    'STEPS', & ! tag name
    'STEPS', & ! fortran variable
    'INTEGER1D', & ! type
    '<NSTP', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_oc_param_definitions(*) = &
    [ &
    gwfoc_budfilerec, &
    gwfoc_budget, &
    gwfoc_fileout, &
    gwfoc_budgetfile, &
    gwfoc_budcsvfilerec, &
    gwfoc_budgetcsv, &
    gwfoc_budgetcsvfile, &
    gwfoc_headfilerec, &
    gwfoc_head, &
    gwfoc_headfile, &
    gwfoc_headprintrec, &
    gwfoc_print_format, &
    gwfoc_formatrecord, &
    gwfoc_columns, &
    gwfoc_width, &
    gwfoc_digits, &
    gwfoc_format, &
    gwfoc_saverecord, &
    gwfoc_save, &
    gwfoc_printrecord, &
    gwfoc_print, &
    gwfoc_rtype, &
    gwfoc_ocsetting, &
    gwfoc_all, &
    gwfoc_first, &
    gwfoc_last, &
    gwfoc_frequency, &
    gwfoc_steps &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_oc_aggregate_definitions(*) = &
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
    gwf_oc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
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

end module GwfOcInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfOcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_oc_param_definitions
  public chf_oc_aggregate_definitions
  public chf_oc_block_definitions
  public ChfOcParamFoundType
  public chf_oc_multi_package
  public chf_oc_subpackages

  type ChfOcParamFoundType
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: qoutfilerec = .false.
    logical :: qoutflow = .false.
    logical :: qoutflowfile = .false.
    logical :: stagefilerec = .false.
    logical :: stage = .false.
    logical :: stagefile = .false.
    logical :: qoutprintrec = .false.
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
  end type ChfOcParamFoundType

  logical :: chf_oc_multi_package = .false.

  character(len=16), parameter :: &
    chf_oc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfoc_budfilerec = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_budget = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_fileout = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_budgetfile = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_budcsvfilerec = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_budgetcsv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_budgetcsvfile = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_qoutfilerec = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'QOUTFLOW_FILERECORD', & ! tag name
    'QOUTFILEREC', & ! fortran variable
    'RECORD QOUTFLOW FILEOUT QOUTFLOWFILE', & ! type
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
    chfoc_qoutflow = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'QOUTFLOW', & ! tag name
    'QOUTFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'qoutflow keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfoc_qoutflowfile = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'QOUTFLOWFILE', & ! tag name
    'QOUTFLOWFILE', & ! fortran variable
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
    chfoc_stagefilerec = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE_FILERECORD', & ! tag name
    'STAGEFILEREC', & ! fortran variable
    'RECORD STAGE FILEOUT STAGEFILE', & ! type
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
    chfoc_stage = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE', & ! tag name
    'STAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'stage keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfoc_stagefile = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'STAGEFILE', & ! tag name
    'STAGEFILE', & ! fortran variable
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
    chfoc_qoutprintrec = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'QOUTFLOWPRINTRECORD', & ! tag name
    'QOUTPRINTREC', & ! fortran variable
    'RECORD QOUTFLOW PRINT_FORMAT FORMATRECORD', & ! type
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
    chfoc_print_format = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_formatrecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_columns = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_width = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_digits = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_format = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_saverecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_save = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_printrecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_print = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_rtype = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_ocsetting = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_all = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_first = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_last = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_frequency = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chfoc_steps = InputParamDefinitionType &
    ( &
    'CHF', & ! component
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
    chf_oc_param_definitions(*) = &
    [ &
    chfoc_budfilerec, &
    chfoc_budget, &
    chfoc_fileout, &
    chfoc_budgetfile, &
    chfoc_budcsvfilerec, &
    chfoc_budgetcsv, &
    chfoc_budgetcsvfile, &
    chfoc_qoutfilerec, &
    chfoc_qoutflow, &
    chfoc_qoutflowfile, &
    chfoc_stagefilerec, &
    chfoc_stage, &
    chfoc_stagefile, &
    chfoc_qoutprintrec, &
    chfoc_print_format, &
    chfoc_formatrecord, &
    chfoc_columns, &
    chfoc_width, &
    chfoc_digits, &
    chfoc_format, &
    chfoc_saverecord, &
    chfoc_save, &
    chfoc_printrecord, &
    chfoc_print, &
    chfoc_rtype, &
    chfoc_ocsetting, &
    chfoc_all, &
    chfoc_first, &
    chfoc_last, &
    chfoc_frequency, &
    chfoc_steps &
    ]

  type(InputParamDefinitionType), parameter :: &
    chf_oc_aggregate_definitions(*) = &
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
    chf_oc_block_definitions(*) = &
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

end module ChfOcInputModule

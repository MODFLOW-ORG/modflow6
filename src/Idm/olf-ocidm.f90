! ** Do Not Modify! MODFLOW 6 system generated file. **
module OlfOcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public olf_oc_param_definitions
  public olf_oc_aggregate_definitions
  public olf_oc_block_definitions
  public OlfOcParamFoundType
  public olf_oc_multi_package
  public olf_oc_subpackages

  type OlfOcParamFoundType
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
  end type OlfOcParamFoundType

  logical :: olf_oc_multi_package = .false.

  character(len=16), parameter :: &
    olf_oc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    olfoc_budfilerec = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_budget = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_fileout = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_budgetfile = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_budcsvfilerec = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_budgetcsv = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_budgetcsvfile = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_qoutfilerec = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_qoutflow = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_qoutflowfile = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_stagefilerec = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_stage = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_stagefile = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_qoutprintrec = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_print_format = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_formatrecord = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_columns = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_width = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_digits = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_format = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_saverecord = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_save = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_printrecord = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_print = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_rtype = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_ocsetting = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_all = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_first = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_last = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_frequency = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olfoc_steps = InputParamDefinitionType &
    ( &
    'OLF', & ! component
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
    olf_oc_param_definitions(*) = &
    [ &
    olfoc_budfilerec, &
    olfoc_budget, &
    olfoc_fileout, &
    olfoc_budgetfile, &
    olfoc_budcsvfilerec, &
    olfoc_budgetcsv, &
    olfoc_budgetcsvfile, &
    olfoc_qoutfilerec, &
    olfoc_qoutflow, &
    olfoc_qoutflowfile, &
    olfoc_stagefilerec, &
    olfoc_stage, &
    olfoc_stagefile, &
    olfoc_qoutprintrec, &
    olfoc_print_format, &
    olfoc_formatrecord, &
    olfoc_columns, &
    olfoc_width, &
    olfoc_digits, &
    olfoc_format, &
    olfoc_saverecord, &
    olfoc_save, &
    olfoc_printrecord, &
    olfoc_print, &
    olfoc_rtype, &
    olfoc_ocsetting, &
    olfoc_all, &
    olfoc_first, &
    olfoc_last, &
    olfoc_frequency, &
    olfoc_steps &
    ]

  type(InputParamDefinitionType), parameter :: &
    olf_oc_aggregate_definitions(*) = &
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
    olf_oc_block_definitions(*) = &
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

end module OlfOcInputModule

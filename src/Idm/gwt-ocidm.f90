! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtOcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_oc_param_definitions
  public gwt_oc_aggregate_definitions
  public gwt_oc_block_definitions
  public GwtOcParamFoundType
  public gwt_oc_multi_package
  public gwt_oc_subpackages

  type GwtOcParamFoundType
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: concfilerec = .false.
    logical :: concentration = .false.
    logical :: concfile = .false.
    logical :: concprintrec = .false.
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
  end type GwtOcParamFoundType

  logical :: gwt_oc_multi_package = .false.

  character(len=16), parameter :: &
    gwt_oc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtoc_budfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_concfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATION_FILERECORD', & ! tag name
    'CONCFILEREC', & ! fortran variable
    'RECORD CONCENTRATION FILEOUT CONCENTRATIONFILE', & ! type
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
    gwtoc_concentration = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'concentration keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtoc_concfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATIONFILE', & ! tag name
    'CONCFILE', & ! fortran variable
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
    gwtoc_concprintrec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'OC', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATIONPRINTRECORD', & ! tag name
    'CONCPRINTREC', & ! fortran variable
    'RECORD CONCENTRATION PRINT_FORMAT FORMATRECORD', & ! type
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
    gwtoc_print_format = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_formatrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_columns = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_width = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_digits = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_format = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_saverecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_save = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_printrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_print = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_rtype = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_ocsetting = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_all = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_first = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_last = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_frequency = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwtoc_steps = InputParamDefinitionType &
    ( &
    'GWT', & ! component
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
    gwt_oc_param_definitions(*) = &
    [ &
    gwtoc_budfilerec, &
    gwtoc_budget, &
    gwtoc_fileout, &
    gwtoc_budgetfile, &
    gwtoc_budcsvfilerec, &
    gwtoc_budgetcsv, &
    gwtoc_budgetcsvfile, &
    gwtoc_concfilerec, &
    gwtoc_concentration, &
    gwtoc_concfile, &
    gwtoc_concprintrec, &
    gwtoc_print_format, &
    gwtoc_formatrecord, &
    gwtoc_columns, &
    gwtoc_width, &
    gwtoc_digits, &
    gwtoc_format, &
    gwtoc_saverecord, &
    gwtoc_save, &
    gwtoc_printrecord, &
    gwtoc_print, &
    gwtoc_rtype, &
    gwtoc_ocsetting, &
    gwtoc_all, &
    gwtoc_first, &
    gwtoc_last, &
    gwtoc_frequency, &
    gwtoc_steps &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwt_oc_aggregate_definitions(*) = &
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
    gwt_oc_block_definitions(*) = &
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

end module GwtOcInputModule

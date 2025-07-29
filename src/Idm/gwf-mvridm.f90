! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfMvrInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_mvr_param_definitions
  public gwf_mvr_aggregate_definitions
  public gwf_mvr_block_definitions
  public GwfMvrParamFoundType
  public gwf_mvr_multi_package
  public gwf_mvr_subpackages

  type GwfMvrParamFoundType
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: modelnames = .false.
    logical :: budgetfr = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budgetcsvfr = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: maxbound = .false.
    logical :: npackages = .false.
    logical :: mname = .false.
    logical :: pname = .false.
    logical :: mname1 = .false.
    logical :: pname1 = .false.
    logical :: id1 = .false.
    logical :: mname2 = .false.
    logical :: pname2 = .false.
    logical :: id2 = .false.
    logical :: mvrtype = .false.
    logical :: value = .false.
  end type GwfMvrParamFoundType

  logical :: gwf_mvr_multi_package = .false.

  character(len=16), parameter :: &
    gwf_mvr_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_print_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_modelnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'MODELNAMES', & ! tag name
    'MODELNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'precede all package names with model names', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budgetfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET_FILERECORD', & ! tag name
    'BUDGETFR', & ! fortran variable
    'RECORD BUDGET FILEOUT BUDGETFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGET', & ! tag name
    'BUDGET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETFILE', & ! tag name
    'BUDGETFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budgetcsvfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV_FILERECORD', & ! tag name
    'BUDGETCSVFR', & ! fortran variable
    'RECORD BUDGETCSV FILEOUT BUDGETCSVFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSV', & ! tag name
    'BUDGETCSV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'OPTIONS', & ! block
    'BUDGETCSVFILE', & ! tag name
    'BUDGETCSVFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXMVR', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of movers', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_npackages = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXPACKAGES', & ! tag name
    'NPACKAGES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of packages to be used with the mover', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_mname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PACKAGES', & ! block
    'MNAME', & ! tag name
    'MNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_pname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PACKAGES', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_mname1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'MNAME1', & ! tag name
    'MNAME1', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_pname1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'PNAME1', & ! tag name
    'PNAME1', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'provider package name', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_id1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'ID1', & ! tag name
    'ID1', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'provider reach', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_mname2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'MNAME2', & ! tag name
    'MNAME2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_pname2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'PNAME2', & ! tag name
    'PNAME2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'receiver package name', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_id2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'ID2', & ! tag name
    'ID2', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'receiver reach', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_mvrtype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'MVRTYPE', & ! tag name
    'MVRTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'mover type', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_value = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'VALUE', & ! tag name
    'VALUE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'mover value', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_mvr_param_definitions(*) = &
    [ &
    gwfmvr_print_input, &
    gwfmvr_print_flows, &
    gwfmvr_modelnames, &
    gwfmvr_budgetfr, &
    gwfmvr_budget, &
    gwfmvr_fileout, &
    gwfmvr_budgetfile, &
    gwfmvr_budgetcsvfr, &
    gwfmvr_budgetcsv, &
    gwfmvr_budgetcsvfile, &
    gwfmvr_maxbound, &
    gwfmvr_npackages, &
    gwfmvr_mname, &
    gwfmvr_pname, &
    gwfmvr_mname1, &
    gwfmvr_pname1, &
    gwfmvr_id1, &
    gwfmvr_mname2, &
    gwfmvr_pname2, &
    gwfmvr_id2, &
    gwfmvr_mvrtype, &
    gwfmvr_value &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_packages = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PACKAGES', & ! block
    'PACKAGES', & ! tag name
    'PACKAGES', & ! fortran variable
    'RECARRAY MNAME PNAME', & ! type
    'NPACKAGES', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfmvr_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'MVR', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY MNAME1 PNAME1 ID1 MNAME2 PNAME2 ID2 MVRTYPE VALUE', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_mvr_aggregate_definitions(*) = &
    [ &
    gwfmvr_packages, &
    gwfmvr_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_mvr_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIMENSIONS', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGES', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfMvrInputModule

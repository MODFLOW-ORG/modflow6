! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtLktInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_lkt_param_definitions
  public gwt_lkt_aggregate_definitions
  public gwt_lkt_block_definitions
  public GwtLktParamFoundType
  public gwt_lkt_multi_package
  public gwt_lkt_subpackages

  type GwtLktParamFoundType
    logical :: flow_pkg_name = .false.
    logical :: auxiliary = .false.
    logical :: fp_aux_name = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprconc = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: conc_filerec = .false.
    logical :: concentration = .false.
    logical :: concfile = .false.
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: dev_nonexpanding = .false.
    logical :: packagedata_ifno = .false.
    logical :: strt = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: concentration_in = .false.
    logical :: rainfall_in = .false.
    logical :: evaporation_in = .false.
    logical :: runoff_in = .false.
    logical :: ext_inflow_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwtLktParamFoundType

  logical :: gwt_lkt_multi_package = .true.

  character(len=16), parameter :: &
    gwt_lkt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_PACKAGE_NAME', & ! tag name
    'FLOW_PKG_NAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of corresponding flow package', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_PACKAGE_AUXILIARY_NAME', & ! tag name
    'FP_AUX_NAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of concentration auxiliary variable '// &
    'in flow package', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_boundnames = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
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
    gwtlkt_iprpak = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
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
    gwtlkt_iprconc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_CONCENTRATION', & ! tag name
    'IPRCONC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated stages to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
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
    gwtlkt_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save lake flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_conc_filerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATION_FILERECORD', & ! tag name
    'CONC_FILEREC', & ! fortran variable
    'RECORD CONCENTRATION FILEOUT CONCFILE', & ! type
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
    gwtlkt_concentration = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION', & ! fortran variable
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
    gwtlkt_concfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'CONCFILE', & ! tag name
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
    gwtlkt_budfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
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
    gwtlkt_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'TS_FILERECORD', & ! tag name
    'TS_FILERECORD', & ! fortran variable
    'RECORD TS6 FILEIN TS6_FILENAME', & ! type
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
    gwtlkt_ts6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'TS6', & ! tag name
    'TS6', & ! fortran variable
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
    gwtlkt_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
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
    gwtlkt_ts6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'TS6_FILENAME', & ! tag name
    'TS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of time series information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
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
    gwtlkt_obs6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_obs6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NONEXPANDING_MATRIX', & ! tag name
    'DEV_NONEXPANDING', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not add rows to the solution matrix', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_strt = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting lake concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_aux = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'AUX', & ! tag name
    'AUX', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_boundname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'lake number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_status = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake concentration status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_concentration_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_rainfall_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'RAINFALL', & ! tag name
    'RAINFALL_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'rainfall concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_evaporation_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'EVAPORATION', & ! tag name
    'EVAPORATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'evaporation concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_runoff_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'RUNOFF', & ! tag name
    'RUNOFF_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'runoff concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_ext_inflow_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'EXT-INFLOW', & ! tag name
    'EXT_INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'ext-inflow concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARYRECORD', & ! tag name
    'AUXILIARYRECORD', & ! fortran variable
    'RECORD AUXILIARY AUXNAME AUXVAL', & ! type
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
    gwtlkt_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'AUXILIARY', & ! tag name
    'PERIOD_AUXILIARY', & ! fortran variable
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
    gwtlkt_auxname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'AUXNAME', & ! tag name
    'AUXNAME', & ! fortran variable
    'STRING', & ! type
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
    gwtlkt_auxval = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'AUXVAL', & ! tag name
    'AUXVAL', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'auxiliary variable value', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_lkt_param_definitions(*) = &
    [ &
    gwtlkt_flow_pkg_name, &
    gwtlkt_auxiliary, &
    gwtlkt_fp_aux_name, &
    gwtlkt_boundnames, &
    gwtlkt_iprpak, &
    gwtlkt_iprconc, &
    gwtlkt_iprflow, &
    gwtlkt_ipakcb, &
    gwtlkt_conc_filerec, &
    gwtlkt_concentration, &
    gwtlkt_concfile, &
    gwtlkt_budfilerec, &
    gwtlkt_budget, &
    gwtlkt_fileout, &
    gwtlkt_budgetfile, &
    gwtlkt_budcsvfilerec, &
    gwtlkt_budgetcsv, &
    gwtlkt_budgetcsvfile, &
    gwtlkt_ts_filerecord, &
    gwtlkt_ts6, &
    gwtlkt_filein, &
    gwtlkt_ts6_filename, &
    gwtlkt_obs_filerecord, &
    gwtlkt_obs6, &
    gwtlkt_obs6_filename, &
    gwtlkt_dev_nonexpanding, &
    gwtlkt_packagedata_ifno, &
    gwtlkt_strt, &
    gwtlkt_aux, &
    gwtlkt_boundname, &
    gwtlkt_ifno, &
    gwtlkt_status, &
    gwtlkt_concentration_in, &
    gwtlkt_rainfall_in, &
    gwtlkt_evaporation_in, &
    gwtlkt_runoff_in, &
    gwtlkt_ext_inflow_in, &
    gwtlkt_auxiliaryrecord, &
    gwtlkt_period_auxiliary, &
    gwtlkt_auxname, &
    gwtlkt_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_packagedata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO STRT AUX BOUNDNAME', & ! type
    'NLAKES', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_lakeperioddata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'LAKEPERIODDATA', & ! tag name
    'LAKEPERIODDATA', & ! fortran variable
    'RECARRAY IFNO LAKSETTING', & ! type
    '', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtlkt_laksetting = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'LKT', & ! subcomponent
    'PERIOD', & ! block
    'LAKSETTING', & ! tag name
    'LAKSETTING', & ! fortran variable
    'KEYSTRING STATUS CONCENTRATION RAINFALL EVAPORATION RUNOFF '// &
    'EXT-INFLOW AUXILIARYRECORD', & ! type
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
    gwt_lkt_aggregate_definitions(*) = &
    [ &
    gwtlkt_packagedata, &
    gwtlkt_lakeperioddata, &
    gwtlkt_laksetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_lkt_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGEDATA', & ! blockname
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

end module GwtLktInputModule

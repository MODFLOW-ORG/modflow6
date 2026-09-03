! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtSftInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_sft_param_definitions
  public gwt_sft_aggregate_definitions
  public gwt_sft_block_definitions
  public GwtSftParamFoundType
  public gwt_sft_multi_package
  public gwt_sft_subpackages

  type GwtSftParamFoundType
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
    logical :: inflow_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwtSftParamFoundType

  logical :: gwt_sft_multi_package = .true.

  character(len=16), parameter :: &
    gwt_sft_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtsft_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_boundnames = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_iprpak = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_iprconc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_CONCENTRATION', & ! tag name
    'IPRCONC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated concentrations to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save reach flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_conc_filerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_concentration = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_concfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_ts6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_ts6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_obs6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_obs6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'reach number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_strt = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting reach concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_aux = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_boundname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'reach number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_status = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'reach concentration status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_concentration_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'reach concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_rainfall_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_evaporation_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_runoff_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_inflow_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'INFLOW', & ! tag name
    'INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'inflow concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_auxname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwtsft_auxval = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
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
    gwt_sft_param_definitions(*) = &
    [ &
    gwtsft_flow_pkg_name, &
    gwtsft_auxiliary, &
    gwtsft_fp_aux_name, &
    gwtsft_boundnames, &
    gwtsft_iprpak, &
    gwtsft_iprconc, &
    gwtsft_iprflow, &
    gwtsft_ipakcb, &
    gwtsft_conc_filerec, &
    gwtsft_concentration, &
    gwtsft_concfile, &
    gwtsft_budfilerec, &
    gwtsft_budget, &
    gwtsft_fileout, &
    gwtsft_budgetfile, &
    gwtsft_budcsvfilerec, &
    gwtsft_budgetcsv, &
    gwtsft_budgetcsvfile, &
    gwtsft_ts_filerecord, &
    gwtsft_ts6, &
    gwtsft_filein, &
    gwtsft_ts6_filename, &
    gwtsft_obs_filerecord, &
    gwtsft_obs6, &
    gwtsft_obs6_filename, &
    gwtsft_dev_nonexpanding, &
    gwtsft_packagedata_ifno, &
    gwtsft_strt, &
    gwtsft_aux, &
    gwtsft_boundname, &
    gwtsft_ifno, &
    gwtsft_status, &
    gwtsft_concentration_in, &
    gwtsft_rainfall_in, &
    gwtsft_evaporation_in, &
    gwtsft_runoff_in, &
    gwtsft_inflow_in, &
    gwtsft_auxiliaryrecord, &
    gwtsft_period_auxiliary, &
    gwtsft_auxname, &
    gwtsft_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtsft_packagedata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO STRT AUX BOUNDNAME', & ! type
    'NREACHES', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtsft_reachperioddata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'REACHPERIODDATA', & ! tag name
    'REACHPERIODDATA', & ! fortran variable
    'RECARRAY IFNO REACHSETTING', & ! type
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
    gwtsft_reachsetting = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'SFT', & ! subcomponent
    'PERIOD', & ! block
    'REACHSETTING', & ! tag name
    'REACHSETTING', & ! fortran variable
    'KEYSTRING STATUS CONCENTRATION RAINFALL EVAPORATION RUNOFF '// &
    'INFLOW AUXILIARYRECORD', & ! type
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
    gwt_sft_aggregate_definitions(*) = &
    [ &
    gwtsft_packagedata, &
    gwtsft_reachperioddata, &
    gwtsft_reachsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_sft_block_definitions(*) = &
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

end module GwtSftInputModule

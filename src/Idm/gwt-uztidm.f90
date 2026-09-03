! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtUztInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_uzt_param_definitions
  public gwt_uzt_aggregate_definitions
  public gwt_uzt_block_definitions
  public GwtUztParamFoundType
  public gwt_uzt_multi_package
  public gwt_uzt_subpackages

  type GwtUztParamFoundType
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
    logical :: infiltration_in = .false.
    logical :: uzet_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwtUztParamFoundType

  logical :: gwt_uzt_multi_package = .true.

  character(len=16), parameter :: &
    gwt_uzt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_boundnames = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_iprpak = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_iprconc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save UZT cell flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_conc_filerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_concentration = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_concfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_ts6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_ts6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_obs6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_obs6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'UZF cell number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_strt = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting UZF cell concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_aux = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_boundname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'UZF cell name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'unsaturated zone flow cell number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_status = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone flow cell concentration status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_concentration_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone flow cell concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_infiltration_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'INFILTRATION', & ! tag name
    'INFILTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'infiltration concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_uzet_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'UZET', & ! tag name
    'UZET_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone et concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_auxname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwtuzt_auxval = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
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
    gwt_uzt_param_definitions(*) = &
    [ &
    gwtuzt_flow_pkg_name, &
    gwtuzt_auxiliary, &
    gwtuzt_fp_aux_name, &
    gwtuzt_boundnames, &
    gwtuzt_iprpak, &
    gwtuzt_iprconc, &
    gwtuzt_iprflow, &
    gwtuzt_ipakcb, &
    gwtuzt_conc_filerec, &
    gwtuzt_concentration, &
    gwtuzt_concfile, &
    gwtuzt_budfilerec, &
    gwtuzt_budget, &
    gwtuzt_fileout, &
    gwtuzt_budgetfile, &
    gwtuzt_budcsvfilerec, &
    gwtuzt_budgetcsv, &
    gwtuzt_budgetcsvfile, &
    gwtuzt_ts_filerecord, &
    gwtuzt_ts6, &
    gwtuzt_filein, &
    gwtuzt_ts6_filename, &
    gwtuzt_obs_filerecord, &
    gwtuzt_obs6, &
    gwtuzt_obs6_filename, &
    gwtuzt_dev_nonexpanding, &
    gwtuzt_packagedata_ifno, &
    gwtuzt_strt, &
    gwtuzt_aux, &
    gwtuzt_boundname, &
    gwtuzt_ifno, &
    gwtuzt_status, &
    gwtuzt_concentration_in, &
    gwtuzt_infiltration_in, &
    gwtuzt_uzet_in, &
    gwtuzt_auxiliaryrecord, &
    gwtuzt_period_auxiliary, &
    gwtuzt_auxname, &
    gwtuzt_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_packagedata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO STRT AUX BOUNDNAME', & ! type
    'NUZFCELLS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtuzt_uztperioddata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'UZTPERIODDATA', & ! tag name
    'UZTPERIODDATA', & ! fortran variable
    'RECARRAY IFNO UZTSETTING', & ! type
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
    gwtuzt_uztsetting = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'UZT', & ! subcomponent
    'PERIOD', & ! block
    'UZTSETTING', & ! tag name
    'UZTSETTING', & ! fortran variable
    'KEYSTRING STATUS CONCENTRATION INFILTRATION UZET '// &
    'AUXILIARYRECORD', & ! type
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
    gwt_uzt_aggregate_definitions(*) = &
    [ &
    gwtuzt_packagedata, &
    gwtuzt_uztperioddata, &
    gwtuzt_uztsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_uzt_block_definitions(*) = &
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

end module GwtUztInputModule

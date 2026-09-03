! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtMwtInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_mwt_param_definitions
  public gwt_mwt_aggregate_definitions
  public gwt_mwt_block_definitions
  public GwtMwtParamFoundType
  public gwt_mwt_multi_package
  public gwt_mwt_subpackages

  type GwtMwtParamFoundType
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
    logical :: rate_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwtMwtParamFoundType

  logical :: gwt_mwt_multi_package = .true.

  character(len=16), parameter :: &
    gwt_mwt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_boundnames = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_iprpak = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_iprconc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_iprflow = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_ipakcb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save well flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_conc_filerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_concentration = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_concfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budget = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budgetfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budgetcsv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_ts6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_ts6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_obs6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_obs6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_strt = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting well concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_aux = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_boundname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_ifno = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'well number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_status = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well concentration status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_concentration_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'CONCENTRATION', & ! tag name
    'CONCENTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_rate_in = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well injection concentration', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_auxname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwtmwt_auxval = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
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
    gwt_mwt_param_definitions(*) = &
    [ &
    gwtmwt_flow_pkg_name, &
    gwtmwt_auxiliary, &
    gwtmwt_fp_aux_name, &
    gwtmwt_boundnames, &
    gwtmwt_iprpak, &
    gwtmwt_iprconc, &
    gwtmwt_iprflow, &
    gwtmwt_ipakcb, &
    gwtmwt_conc_filerec, &
    gwtmwt_concentration, &
    gwtmwt_concfile, &
    gwtmwt_budfilerec, &
    gwtmwt_budget, &
    gwtmwt_fileout, &
    gwtmwt_budgetfile, &
    gwtmwt_budcsvfilerec, &
    gwtmwt_budgetcsv, &
    gwtmwt_budgetcsvfile, &
    gwtmwt_ts_filerecord, &
    gwtmwt_ts6, &
    gwtmwt_filein, &
    gwtmwt_ts6_filename, &
    gwtmwt_obs_filerecord, &
    gwtmwt_obs6, &
    gwtmwt_obs6_filename, &
    gwtmwt_dev_nonexpanding, &
    gwtmwt_packagedata_ifno, &
    gwtmwt_strt, &
    gwtmwt_aux, &
    gwtmwt_boundname, &
    gwtmwt_ifno, &
    gwtmwt_status, &
    gwtmwt_concentration_in, &
    gwtmwt_rate_in, &
    gwtmwt_auxiliaryrecord, &
    gwtmwt_period_auxiliary, &
    gwtmwt_auxname, &
    gwtmwt_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_packagedata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO STRT AUX BOUNDNAME', & ! type
    'NMAWWELLS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtmwt_mwtperioddata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'MWTPERIODDATA', & ! tag name
    'MWTPERIODDATA', & ! fortran variable
    'RECARRAY IFNO MWTSETTING', & ! type
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
    gwtmwt_mwtsetting = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'MWT', & ! subcomponent
    'PERIOD', & ! block
    'MWTSETTING', & ! tag name
    'MWTSETTING', & ! fortran variable
    'KEYSTRING STATUS CONCENTRATION RATE AUXILIARYRECORD', & ! type
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
    gwt_mwt_aggregate_definitions(*) = &
    [ &
    gwtmwt_packagedata, &
    gwtmwt_mwtperioddata, &
    gwtmwt_mwtsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_mwt_block_definitions(*) = &
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

end module GwtMwtInputModule

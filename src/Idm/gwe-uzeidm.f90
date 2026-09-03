! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweUzeInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_uze_param_definitions
  public gwe_uze_aggregate_definitions
  public gwe_uze_block_definitions
  public GweUzeParamFoundType
  public gwe_uze_multi_package
  public gwe_uze_subpackages

  type GweUzeParamFoundType
    logical :: flow_pkg_name = .false.
    logical :: auxiliary = .false.
    logical :: fp_aux_name = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprconc = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: temp_filerec = .false.
    logical :: temperature = .false.
    logical :: tempfile = .false.
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
    logical :: temperature_in = .false.
    logical :: infiltration_in = .false.
    logical :: uzet_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GweUzeParamFoundType

  logical :: gwe_uze_multi_package = .true.

  character(len=16), parameter :: &
    gwe_uze_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gweuze_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_boundnames = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_iprpak = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_iprconc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_TEMPERATURE', & ! tag name
    'IPRCONC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated temperatures to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_iprflow = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_ipakcb = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save UZE cell flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_temp_filerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'OPTIONS', & ! block
    'TEMPERATURE_FILERECORD', & ! tag name
    'TEMP_FILEREC', & ! fortran variable
    'RECORD TEMPERATURE FILEOUT TEMPFILE', & ! type
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
    gweuze_temperature = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'OPTIONS', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE', & ! fortran variable
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
    gweuze_tempfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'OPTIONS', & ! block
    'TEMPFILE', & ! tag name
    'TEMPFILE', & ! fortran variable
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
    gweuze_budfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_budget = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_budgetfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_budgetcsv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_ts6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_ts6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_obs6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_obs6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'UZFNO', & ! tag name
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
    gweuze_strt = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting UZF cell temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_aux = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_boundname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'UZFNO', & ! tag name
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
    gweuze_status = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone flow cell temperature status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_temperature_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone flow cell temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_infiltration_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'INFILTRATION', & ! tag name
    'INFILTRATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'infiltration temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_uzet_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'UZET', & ! tag name
    'UZET_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'unsaturated zone et temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gweuze_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_auxname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gweuze_auxval = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
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
    gwe_uze_param_definitions(*) = &
    [ &
    gweuze_flow_pkg_name, &
    gweuze_auxiliary, &
    gweuze_fp_aux_name, &
    gweuze_boundnames, &
    gweuze_iprpak, &
    gweuze_iprconc, &
    gweuze_iprflow, &
    gweuze_ipakcb, &
    gweuze_temp_filerec, &
    gweuze_temperature, &
    gweuze_tempfile, &
    gweuze_budfilerec, &
    gweuze_budget, &
    gweuze_fileout, &
    gweuze_budgetfile, &
    gweuze_budcsvfilerec, &
    gweuze_budgetcsv, &
    gweuze_budgetcsvfile, &
    gweuze_ts_filerecord, &
    gweuze_ts6, &
    gweuze_filein, &
    gweuze_ts6_filename, &
    gweuze_obs_filerecord, &
    gweuze_obs6, &
    gweuze_obs6_filename, &
    gweuze_dev_nonexpanding, &
    gweuze_packagedata_ifno, &
    gweuze_strt, &
    gweuze_aux, &
    gweuze_boundname, &
    gweuze_ifno, &
    gweuze_status, &
    gweuze_temperature_in, &
    gweuze_infiltration_in, &
    gweuze_uzet_in, &
    gweuze_auxiliaryrecord, &
    gweuze_period_auxiliary, &
    gweuze_auxname, &
    gweuze_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gweuze_packagedata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY UZFNO STRT AUX BOUNDNAME', & ! type
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
    gweuze_uzeperioddata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'UZEPERIODDATA', & ! tag name
    'UZEPERIODDATA', & ! fortran variable
    'RECARRAY UZFNO UZESETTING', & ! type
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
    gweuze_uzesetting = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'UZE', & ! subcomponent
    'PERIOD', & ! block
    'UZESETTING', & ! tag name
    'UZESETTING', & ! fortran variable
    'KEYSTRING STATUS TEMPERATURE INFILTRATION UZET '// &
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
    gwe_uze_aggregate_definitions(*) = &
    [ &
    gweuze_packagedata, &
    gweuze_uzeperioddata, &
    gweuze_uzesetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_uze_block_definitions(*) = &
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

end module GweUzeInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweMweInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_mwe_param_definitions
  public gwe_mwe_aggregate_definitions
  public gwe_mwe_block_definitions
  public GweMweParamFoundType
  public gwe_mwe_multi_package
  public gwe_mwe_subpackages

  type GweMweParamFoundType
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
    logical :: ktf = .false.
    logical :: fthk = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: temperature_in = .false.
    logical :: rate_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GweMweParamFoundType

  logical :: gwe_mwe_multi_package = .true.

  character(len=16), parameter :: &
    gwe_mwe_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwemwe_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_PACKAGE_AUXILIARY_NAME', & ! tag name
    'FP_AUX_NAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of temperature auxiliary variable in '// &
    'flow package', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_boundnames = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_iprpak = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_iprconc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_iprflow = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_ipakcb = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_temp_filerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_temperature = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_tempfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budget = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budgetfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budgetcsv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_ts6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_ts6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_obs6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_obs6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'MAWNO', & ! tag name
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
    gwemwe_strt = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting well temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_ktf = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'KTF', & ! tag name
    'KTF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'thermal conductivity of the feature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_fthk = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'FTHK', & ! tag name
    'FTHK', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'thickness of the well feature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_aux = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_boundname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'MAWNO', & ! tag name
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
    gwemwe_status = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well temperature status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_temperature_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_rate_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'RATE', & ! tag name
    'RATE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'well injection temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwemwe_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_auxname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwemwe_auxval = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
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
    gwe_mwe_param_definitions(*) = &
    [ &
    gwemwe_flow_pkg_name, &
    gwemwe_auxiliary, &
    gwemwe_fp_aux_name, &
    gwemwe_boundnames, &
    gwemwe_iprpak, &
    gwemwe_iprconc, &
    gwemwe_iprflow, &
    gwemwe_ipakcb, &
    gwemwe_temp_filerec, &
    gwemwe_temperature, &
    gwemwe_tempfile, &
    gwemwe_budfilerec, &
    gwemwe_budget, &
    gwemwe_fileout, &
    gwemwe_budgetfile, &
    gwemwe_budcsvfilerec, &
    gwemwe_budgetcsv, &
    gwemwe_budgetcsvfile, &
    gwemwe_ts_filerecord, &
    gwemwe_ts6, &
    gwemwe_filein, &
    gwemwe_ts6_filename, &
    gwemwe_obs_filerecord, &
    gwemwe_obs6, &
    gwemwe_obs6_filename, &
    gwemwe_dev_nonexpanding, &
    gwemwe_packagedata_ifno, &
    gwemwe_strt, &
    gwemwe_ktf, &
    gwemwe_fthk, &
    gwemwe_aux, &
    gwemwe_boundname, &
    gwemwe_ifno, &
    gwemwe_status, &
    gwemwe_temperature_in, &
    gwemwe_rate_in, &
    gwemwe_auxiliaryrecord, &
    gwemwe_period_auxiliary, &
    gwemwe_auxname, &
    gwemwe_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwemwe_packagedata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY MAWNO STRT KTF FTHK AUX BOUNDNAME', & ! type
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
    gwemwe_mweperioddata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'MWEPERIODDATA', & ! tag name
    'MWEPERIODDATA', & ! fortran variable
    'RECARRAY MAWNO MWESETTING', & ! type
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
    gwemwe_mwesetting = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'MWE', & ! subcomponent
    'PERIOD', & ! block
    'MWESETTING', & ! tag name
    'MWESETTING', & ! fortran variable
    'KEYSTRING STATUS TEMPERATURE RATE AUXILIARYRECORD', & ! type
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
    gwe_mwe_aggregate_definitions(*) = &
    [ &
    gwemwe_packagedata, &
    gwemwe_mweperioddata, &
    gwemwe_mwesetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_mwe_block_definitions(*) = &
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

end module GweMweInputModule

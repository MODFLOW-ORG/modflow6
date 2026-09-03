! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweLkeInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_lke_param_definitions
  public gwe_lke_aggregate_definitions
  public gwe_lke_block_definitions
  public GweLkeParamFoundType
  public gwe_lke_multi_package
  public gwe_lke_subpackages

  type GweLkeParamFoundType
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
    logical :: rbthcnd = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: temperature_in = .false.
    logical :: rainfall_in = .false.
    logical :: evaporation_in = .false.
    logical :: runoff_in = .false.
    logical :: ext_inflow_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GweLkeParamFoundType

  logical :: gwe_lke_multi_package = .true.

  character(len=16), parameter :: &
    gwe_lke_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwelke_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_boundnames = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_iprpak = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_iprconc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_iprflow = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_ipakcb = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_temp_filerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_temperature = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_tempfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budget = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budgetfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budgetcsv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_ts6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_ts6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_obs6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_obs6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'LAKENO', & ! tag name
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
    gwelke_strt = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting lake temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_ktf = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'KTF', & ! tag name
    'KTF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'boundary thermal conductivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_rbthcnd = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RBTHCND', & ! tag name
    'RBTHCND', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'streambed thickness', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_aux = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_boundname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'LAKENO', & ! tag name
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
    gwelke_status = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake temperature status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_temperature_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'lake temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_rainfall_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'RAINFALL', & ! tag name
    'RAINFALL_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'rainfall temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_evaporation_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'EVAPORATION', & ! tag name
    'EVAPORATION_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'evaporation temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_runoff_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'RUNOFF', & ! tag name
    'RUNOFF_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'runoff temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_ext_inflow_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'EXT-INFLOW', & ! tag name
    'EXT_INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'ext-inflow temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwelke_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_auxname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwelke_auxval = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
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
    gwe_lke_param_definitions(*) = &
    [ &
    gwelke_flow_pkg_name, &
    gwelke_auxiliary, &
    gwelke_fp_aux_name, &
    gwelke_boundnames, &
    gwelke_iprpak, &
    gwelke_iprconc, &
    gwelke_iprflow, &
    gwelke_ipakcb, &
    gwelke_temp_filerec, &
    gwelke_temperature, &
    gwelke_tempfile, &
    gwelke_budfilerec, &
    gwelke_budget, &
    gwelke_fileout, &
    gwelke_budgetfile, &
    gwelke_budcsvfilerec, &
    gwelke_budgetcsv, &
    gwelke_budgetcsvfile, &
    gwelke_ts_filerecord, &
    gwelke_ts6, &
    gwelke_filein, &
    gwelke_ts6_filename, &
    gwelke_obs_filerecord, &
    gwelke_obs6, &
    gwelke_obs6_filename, &
    gwelke_dev_nonexpanding, &
    gwelke_packagedata_ifno, &
    gwelke_strt, &
    gwelke_ktf, &
    gwelke_rbthcnd, &
    gwelke_aux, &
    gwelke_boundname, &
    gwelke_ifno, &
    gwelke_status, &
    gwelke_temperature_in, &
    gwelke_rainfall_in, &
    gwelke_evaporation_in, &
    gwelke_runoff_in, &
    gwelke_ext_inflow_in, &
    gwelke_auxiliaryrecord, &
    gwelke_period_auxiliary, &
    gwelke_auxname, &
    gwelke_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwelke_packagedata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY LAKENO STRT KTF RBTHCND AUX BOUNDNAME', & ! type
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
    gwelke_lakeperioddata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'LAKEPERIODDATA', & ! tag name
    'LAKEPERIODDATA', & ! fortran variable
    'RECARRAY LAKENO LAKSETTING', & ! type
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
    gwelke_laksetting = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'LKE', & ! subcomponent
    'PERIOD', & ! block
    'LAKSETTING', & ! tag name
    'LAKSETTING', & ! fortran variable
    'KEYSTRING STATUS TEMPERATURE RAINFALL EVAPORATION RUNOFF '// &
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
    gwe_lke_aggregate_definitions(*) = &
    [ &
    gwelke_packagedata, &
    gwelke_lakeperioddata, &
    gwelke_laksetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_lke_block_definitions(*) = &
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

end module GweLkeInputModule

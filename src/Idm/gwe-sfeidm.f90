! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweSfeInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_sfe_param_definitions
  public gwe_sfe_aggregate_definitions
  public gwe_sfe_block_definitions
  public GweSfeParamFoundType
  public gwe_sfe_multi_package
  public gwe_sfe_subpackages

  type GweSfeParamFoundType
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
    logical :: inflow_in = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GweSfeParamFoundType

  logical :: gwe_sfe_multi_package = .true.

  character(len=16), parameter :: &
    gwe_sfe_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwesfe_flow_pkg_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_fp_aux_name = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_boundnames = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_iprpak = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_iprconc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_TEMPERATURE', & ! tag name
    'IPRCONC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated temperature to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_iprflow = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_ipakcb = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_temp_filerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_temperature = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'OPTIONS', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'temperature keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_tempfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budget = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budgetfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budgetcsv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_ts6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_ts6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_obs6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_obs6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_dev_nonexpanding = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RNO', & ! tag name
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
    gwesfe_strt = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'STRT', & ! tag name
    'STRT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'starting reach temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_ktf = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_rbthcnd = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_aux = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_boundname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_ifno = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'RNO', & ! tag name
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
    gwesfe_status = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'reach temperature status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_temperature_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'TEMPERATURE', & ! tag name
    'TEMPERATURE_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'reach temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_rainfall_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_evaporation_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_runoff_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_inflow_in = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'INFLOW', & ! tag name
    'INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'inflow temperature', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwesfe_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_auxname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwesfe_auxval = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
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
    gwe_sfe_param_definitions(*) = &
    [ &
    gwesfe_flow_pkg_name, &
    gwesfe_auxiliary, &
    gwesfe_fp_aux_name, &
    gwesfe_boundnames, &
    gwesfe_iprpak, &
    gwesfe_iprconc, &
    gwesfe_iprflow, &
    gwesfe_ipakcb, &
    gwesfe_temp_filerec, &
    gwesfe_temperature, &
    gwesfe_tempfile, &
    gwesfe_budfilerec, &
    gwesfe_budget, &
    gwesfe_fileout, &
    gwesfe_budgetfile, &
    gwesfe_budcsvfilerec, &
    gwesfe_budgetcsv, &
    gwesfe_budgetcsvfile, &
    gwesfe_ts_filerecord, &
    gwesfe_ts6, &
    gwesfe_filein, &
    gwesfe_ts6_filename, &
    gwesfe_obs_filerecord, &
    gwesfe_obs6, &
    gwesfe_obs6_filename, &
    gwesfe_dev_nonexpanding, &
    gwesfe_packagedata_ifno, &
    gwesfe_strt, &
    gwesfe_ktf, &
    gwesfe_rbthcnd, &
    gwesfe_aux, &
    gwesfe_boundname, &
    gwesfe_ifno, &
    gwesfe_status, &
    gwesfe_temperature_in, &
    gwesfe_rainfall_in, &
    gwesfe_evaporation_in, &
    gwesfe_runoff_in, &
    gwesfe_inflow_in, &
    gwesfe_auxiliaryrecord, &
    gwesfe_period_auxiliary, &
    gwesfe_auxname, &
    gwesfe_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwesfe_packagedata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY RNO STRT KTF RBTHCND AUX BOUNDNAME', & ! type
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
    gwesfe_reachperioddata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'REACHPERIODDATA', & ! tag name
    'REACHPERIODDATA', & ! fortran variable
    'RECARRAY RNO REACHSETTING', & ! type
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
    gwesfe_reachsetting = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'SFE', & ! subcomponent
    'PERIOD', & ! block
    'REACHSETTING', & ! tag name
    'REACHSETTING', & ! fortran variable
    'KEYSTRING STATUS TEMPERATURE RAINFALL EVAPORATION RUNOFF '// &
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
    gwe_sfe_aggregate_definitions(*) = &
    [ &
    gwesfe_packagedata, &
    gwesfe_reachperioddata, &
    gwesfe_reachsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_sfe_block_definitions(*) = &
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

end module GweSfeInputModule

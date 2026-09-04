! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfSfrInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_sfr_param_definitions
  public gwf_sfr_aggregate_definitions
  public gwf_sfr_block_definitions
  public GwfSfrParamFoundType
  public gwf_sfr_multi_package
  public gwf_sfr_subpackages

  type GwfSfrParamFoundType
    logical :: storage = .false.
    logical :: ats_courant = .false.
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: print_stage = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: stage_filerecord = .false.
    logical :: options_stage = .false.
    logical :: stagefile = .false.
    logical :: budfilerec = .false.
    logical :: budget = .false.
    logical :: fileout = .false.
    logical :: budgetfile = .false.
    logical :: budcsvfilerec = .false.
    logical :: budgetcsv = .false.
    logical :: budgetcsvfile = .false.
    logical :: pkgconvfilerec = .false.
    logical :: pkgconv = .false.
    logical :: pkgconvfname = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: options_filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: maxsfrpicard = .false.
    logical :: maxsfrit = .false.
    logical :: dmaxchg = .false.
    logical :: unit_conversion = .false.
    logical :: lengthconv = .false.
    logical :: time_conversion = .false.
    logical :: storage_weight = .false.
    logical :: icheck = .false.
    logical :: iconvchk = .false.
    logical :: nreaches = .false.
    logical :: packagedata_ifno = .false.
    logical :: cellid = .false.
    logical :: rlen = .false.
    logical :: rwid = .false.
    logical :: rgrd = .false.
    logical :: rtp = .false.
    logical :: rbth = .false.
    logical :: rhk = .false.
    logical :: man = .false.
    logical :: ncon = .false.
    logical :: ustrf = .false.
    logical :: ndv = .false.
    logical :: aux = .false.
    logical :: boundname = .false.
    logical :: crosssect_ifno = .false.
    logical :: xs_tab6 = .false.
    logical :: xs_filein = .false.
    logical :: xs_tab6_filename = .false.
    logical :: conndata_ifno = .false.
    logical :: ic = .false.
    logical :: diversions_ifno = .false.
    logical :: diversions_idv = .false.
    logical :: iconr = .false.
    logical :: cprior = .false.
    logical :: initstage_ifno = .false.
    logical :: initialstage = .false.
    logical :: ifno = .false.
    logical :: status = .false.
    logical :: bedk = .false.
    logical :: manning = .false.
    logical :: stage_in = .false.
    logical :: inflow_in = .false.
    logical :: rainfall_in = .false.
    logical :: evaporation_in = .false.
    logical :: runoff_in = .false.
    logical :: diversionrecord = .false.
    logical :: diversion = .false.
    logical :: idv = .false.
    logical :: divflow = .false.
    logical :: upstream_frac = .false.
    logical :: xsectionrecord = .false.
    logical :: cross_section = .false.
    logical :: tab6 = .false.
    logical :: filein = .false.
    logical :: tab6_filename = .false.
    logical :: auxiliaryrecord = .false.
    logical :: period_auxiliary = .false.
    logical :: auxname = .false.
    logical :: auxval = .false.
  end type GwfSfrParamFoundType

  logical :: gwf_sfr_multi_package = .true.

  character(len=16), parameter :: &
    gwf_sfr_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_storage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'STORAGE', & ! tag name
    'STORAGE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate reach storage', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_ats_courant = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'ATS_COURANT', & ! tag name
    'ATS_COURANT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'target Courant number for adaptive time stepping', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_print_stage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_STAGE', & ! tag name
    'PRINT_STAGE', & ! fortran variable
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
    gwfsfr_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save stream reach flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_stage_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE_FILERECORD', & ! tag name
    'STAGE_FILERECORD', & ! fortran variable
    'RECORD STAGE FILEOUT STAGEFILE', & ! type
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
    gwfsfr_options_stage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'STAGE', & ! tag name
    'OPTIONS_STAGE', & ! fortran variable
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
    gwfsfr_stagefile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'STAGEFILE', & ! tag name
    'STAGEFILE', & ! fortran variable
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
    gwfsfr_budfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_pkgconvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILERECORD', & ! tag name
    'PKGCONVFILEREC', & ! fortran variable
    'RECORD PACKAGE_CONVERGENCE FILEOUT '// &
    'PACKAGE_CONVERGENCE_FILENAME', & ! type
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
    gwfsfr_pkgconv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE', & ! tag name
    'PKGCONV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'package_convergence keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_pkgconvfname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILENAME', & ! tag name
    'PKGCONVFNAME', & ! fortran variable
    'STRING', & ! type
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
    gwfsfr_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_options_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'OPTIONS_FILEIN', & ! fortran variable
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
    gwfsfr_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'MOVER', & ! tag name
    'MOVER', & ! fortran variable
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
    gwfsfr_maxsfrpicard = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'MAXIMUM_PICARD_ITERATIONS', & ! tag name
    'MAXSFRPICARD', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'SFR picard iterations', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_maxsfrit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'MAXIMUM_ITERATIONS', & ! tag name
    'MAXSFRIT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'SFR Newton-Raphson iterations', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_dmaxchg = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'MAXIMUM_DEPTH_CHANGE', & ! tag name
    'DMAXCHG', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'depth closure tolerance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_unit_conversion = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'UNIT_CONVERSION', & ! tag name
    'UNIT_CONVERSION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_lengthconv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_CONVERSION', & ! tag name
    'LENGTHCONV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'length conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_time_conversion = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'TIME_CONVERSION', & ! tag name
    'TIME_CONVERSION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'time conversion factor', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_storage_weight = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_STORAGE_WEIGHT', & ! tag name
    'STORAGE_WEIGHT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach storage time weighting', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_icheck = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_CHECK', & ! tag name
    'ICHECK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not check reach geometry and parameters', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_iconvchk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_FINAL_CHECK', & ! tag name
    'ICONVCHK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not perform final convergence check', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_nreaches = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIMENSIONS', & ! block
    'NREACHES', & ! tag name
    'NREACHES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number reaches', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rlen = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RLEN', & ! tag name
    'RLEN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach length', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rwid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RWID', & ! tag name
    'RWID', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rgrd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RGRD', & ! tag name
    'RGRD', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'stream gradient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rtp = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RTP', & ! tag name
    'RTP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach streambed top elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rbth = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RBTH', & ! tag name
    'RBTH', & ! fortran variable
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
    gwfsfr_rhk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RHK', & ! tag name
    'RHK', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach bed hydraulic conductivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_man = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'MAN', & ! tag name
    'MAN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'Mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_ncon = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NCON', & ! tag name
    'NCON', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of connected reaches', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_ustrf = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'USTRF', & ! tag name
    'USTRF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fraction of upstream flow', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_ndv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NDV', & ! tag name
    'NDV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of downstream reaches', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_aux = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_crosssect_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CROSSSECTIONS', & ! block
    'IFNO', & ! tag name
    'CROSSSECT_IFNO', & ! fortran variable
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
    gwfsfr_xs_tab6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CROSSSECTIONS', & ! block
    'TAB6', & ! tag name
    'XS_TAB6', & ! fortran variable
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
    gwfsfr_xs_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CROSSSECTIONS', & ! block
    'FILEIN', & ! tag name
    'XS_FILEIN', & ! fortran variable
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
    gwfsfr_xs_tab6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CROSSSECTIONS', & ! block
    'TAB6_FILENAME', & ! tag name
    'XS_TAB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'table file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_conndata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IFNO', & ! tag name
    'CONNDATA_IFNO', & ! fortran variable
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
    gwfsfr_ic = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IC', & ! tag name
    'IC', & ! fortran variable
    'INTEGER1D', & ! type
    ':', & ! shape
    'connected reach numbers', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_diversions_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIVERSIONS', & ! block
    'IFNO', & ! tag name
    'DIVERSIONS_IFNO', & ! fortran variable
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
    gwfsfr_diversions_idv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIVERSIONS', & ! block
    'IDV', & ! tag name
    'DIVERSIONS_IDV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'downstream diversion number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_iconr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIVERSIONS', & ! block
    'ICONR', & ! tag name
    'ICONR', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'downstream reach number for diversion', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_cprior = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIVERSIONS', & ! block
    'CPRIOR', & ! tag name
    'CPRIOR', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'iprior code', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_initstage_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'INITIALSTAGES', & ! block
    'IFNO', & ! tag name
    'INITSTAGE_IFNO', & ! fortran variable
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
    gwfsfr_initialstage = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'INITIALSTAGES', & ! block
    'INITIALSTAGE', & ! tag name
    'INITIALSTAGE', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial reach stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_status = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'STATUS', & ! tag name
    'STATUS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'reach status', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_bedk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'BEDK', & ! tag name
    'BEDK', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach bed hydraulic conductivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_manning = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'MANNING', & ! tag name
    'MANNING', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reach Mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_stage_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'STAGE', & ! tag name
    'STAGE_IN', & ! fortran variable
    'STRING', & ! type
    'NREACHES', & ! shape
    'reach stage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_inflow_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'INFLOW', & ! tag name
    'INFLOW_IN', & ! fortran variable
    'STRING', & ! type
    'NREACHES', & ! shape
    'inflow rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_rainfall_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'RAINFALL', & ! tag name
    'RAINFALL_IN', & ! fortran variable
    'STRING', & ! type
    'NREACHES', & ! shape
    'rainfall rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_evaporation_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'EVAPORATION', & ! tag name
    'EVAPORATION_IN', & ! fortran variable
    'STRING', & ! type
    'NREACHES', & ! shape
    'evaporation rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_runoff_in = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'RUNOFF', & ! tag name
    'RUNOFF_IN', & ! fortran variable
    'STRING', & ! type
    'NREACHES', & ! shape
    'runoff rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_diversionrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'DIVERSIONRECORD', & ! tag name
    'DIVERSIONRECORD', & ! fortran variable
    'RECORD DIVERSION IDV DIVFLOW', & ! type
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
    gwfsfr_diversion = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'DIVERSION', & ! tag name
    'DIVERSION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'diversion keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_idv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'IDV', & ! tag name
    'IDV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'diversion number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_divflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'DIVFLOW', & ! tag name
    'DIVFLOW', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'volumetric diversion flow rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_upstream_frac = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'UPSTREAM_FRACTION', & ! tag name
    'UPSTREAM_FRAC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fraction of upstream flow', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_xsectionrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'CROSS_SECTIONRECORD', & ! tag name
    'XSECTIONRECORD', & ! fortran variable
    'RECORD CROSS_SECTION TAB6 FILEIN TAB6_FILENAME', & ! type
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
    gwfsfr_cross_section = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'CROSS_SECTION', & ! tag name
    'CROSS_SECTION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'cross_section keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_tab6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'TAB6', & ! tag name
    'TAB6', & ! fortran variable
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
    gwfsfr_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
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
    gwfsfr_tab6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'TAB6_FILENAME', & ! tag name
    'TAB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'table file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_auxiliaryrecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_period_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_auxname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwfsfr_auxval = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
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
    gwf_sfr_param_definitions(*) = &
    [ &
    gwfsfr_storage, &
    gwfsfr_ats_courant, &
    gwfsfr_auxiliary, &
    gwfsfr_boundnames, &
    gwfsfr_iprpak, &
    gwfsfr_print_stage, &
    gwfsfr_iprflow, &
    gwfsfr_ipakcb, &
    gwfsfr_stage_filerecord, &
    gwfsfr_options_stage, &
    gwfsfr_stagefile, &
    gwfsfr_budfilerec, &
    gwfsfr_budget, &
    gwfsfr_fileout, &
    gwfsfr_budgetfile, &
    gwfsfr_budcsvfilerec, &
    gwfsfr_budgetcsv, &
    gwfsfr_budgetcsvfile, &
    gwfsfr_pkgconvfilerec, &
    gwfsfr_pkgconv, &
    gwfsfr_pkgconvfname, &
    gwfsfr_ts_filerecord, &
    gwfsfr_ts6, &
    gwfsfr_options_filein, &
    gwfsfr_ts6_filename, &
    gwfsfr_obs_filerecord, &
    gwfsfr_obs6, &
    gwfsfr_obs6_filename, &
    gwfsfr_mover, &
    gwfsfr_maxsfrpicard, &
    gwfsfr_maxsfrit, &
    gwfsfr_dmaxchg, &
    gwfsfr_unit_conversion, &
    gwfsfr_lengthconv, &
    gwfsfr_time_conversion, &
    gwfsfr_storage_weight, &
    gwfsfr_icheck, &
    gwfsfr_iconvchk, &
    gwfsfr_nreaches, &
    gwfsfr_packagedata_ifno, &
    gwfsfr_cellid, &
    gwfsfr_rlen, &
    gwfsfr_rwid, &
    gwfsfr_rgrd, &
    gwfsfr_rtp, &
    gwfsfr_rbth, &
    gwfsfr_rhk, &
    gwfsfr_man, &
    gwfsfr_ncon, &
    gwfsfr_ustrf, &
    gwfsfr_ndv, &
    gwfsfr_aux, &
    gwfsfr_boundname, &
    gwfsfr_crosssect_ifno, &
    gwfsfr_xs_tab6, &
    gwfsfr_xs_filein, &
    gwfsfr_xs_tab6_filename, &
    gwfsfr_conndata_ifno, &
    gwfsfr_ic, &
    gwfsfr_diversions_ifno, &
    gwfsfr_diversions_idv, &
    gwfsfr_iconr, &
    gwfsfr_cprior, &
    gwfsfr_initstage_ifno, &
    gwfsfr_initialstage, &
    gwfsfr_ifno, &
    gwfsfr_status, &
    gwfsfr_bedk, &
    gwfsfr_manning, &
    gwfsfr_stage_in, &
    gwfsfr_inflow_in, &
    gwfsfr_rainfall_in, &
    gwfsfr_evaporation_in, &
    gwfsfr_runoff_in, &
    gwfsfr_diversionrecord, &
    gwfsfr_diversion, &
    gwfsfr_idv, &
    gwfsfr_divflow, &
    gwfsfr_upstream_frac, &
    gwfsfr_xsectionrecord, &
    gwfsfr_cross_section, &
    gwfsfr_tab6, &
    gwfsfr_filein, &
    gwfsfr_tab6_filename, &
    gwfsfr_auxiliaryrecord, &
    gwfsfr_period_auxiliary, &
    gwfsfr_auxname, &
    gwfsfr_auxval &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO CELLID RLEN RWID RGRD RTP RBTH RHK MAN NCON '// &
    'USTRF NDV AUX BOUNDNAME', & ! type
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
    gwfsfr_crosssections = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CROSSSECTIONS', & ! block
    'CROSSSECTIONS', & ! tag name
    'CROSSSECTIONS', & ! fortran variable
    'RECARRAY IFNO TAB6 FILEIN TAB6_FILENAME', & ! type
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
    gwfsfr_connectiondata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CONNECTIONDATA', & ! tag name
    'CONNECTIONDATA', & ! fortran variable
    'RECARRAY IFNO IC', & ! type
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
    gwfsfr_diversions = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'DIVERSIONS', & ! block
    'DIVERSIONS', & ! tag name
    'DIVERSIONS', & ! fortran variable
    'RECARRAY IFNO IDV ICONR CPRIOR', & ! type
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
    gwfsfr_initialstages = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'INITIALSTAGES', & ! block
    'INITIALSTAGES', & ! tag name
    'INITIALSTAGES', & ! fortran variable
    'RECARRAY IFNO INITIALSTAGE', & ! type
    'NREACHES', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsfr_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY IFNO SFRSETTING', & ! type
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
    gwfsfr_sfrsetting = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SFR', & ! subcomponent
    'PERIOD', & ! block
    'SFRSETTING', & ! tag name
    'SFRSETTING', & ! fortran variable
    'KEYSTRING STATUS BEDK MANNING STAGE INFLOW RAINFALL '// &
    'EVAPORATION RUNOFF DIVERSIONRECORD UPSTREAM_FRACTION '// &
    'CROSS_SECTIONRECORD AUXILIARYRECORD', & ! type
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
    gwf_sfr_aggregate_definitions(*) = &
    [ &
    gwfsfr_packagedata, &
    gwfsfr_crosssections, &
    gwfsfr_connectiondata, &
    gwfsfr_diversions, &
    gwfsfr_initialstages, &
    gwfsfr_perioddata, &
    gwfsfr_sfrsetting &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_sfr_block_definitions(*) = &
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
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CROSSSECTIONS', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CONNECTIONDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'DIVERSIONS', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'INITIALSTAGES', & ! blockname
    .false., & ! required
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

end module GwfSfrInputModule

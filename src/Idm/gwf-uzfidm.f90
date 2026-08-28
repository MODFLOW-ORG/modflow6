! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfUzfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_uzf_param_definitions
  public gwf_uzf_aggregate_definitions
  public gwf_uzf_block_definitions
  public GwfUzfParamFoundType
  public gwf_uzf_multi_package
  public gwf_uzf_subpackages

  type GwfUzfParamFoundType
    logical :: auxiliary = .false.
    logical :: auxmultname = .false.
    logical :: boundnames = .false.
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: ipakcb = .false.
    logical :: wc_filerecord = .false.
    logical :: water_content = .false.
    logical :: wcfile = .false.
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
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: mover = .false.
    logical :: iconvchk = .false.
    logical :: simulate_et = .false.
    logical :: linear_gwet = .false.
    logical :: square_gwet = .false.
    logical :: simulate_gwseep = .false.
    logical :: unsat_etwc = .false.
    logical :: unsat_etae = .false.
    logical :: maxbound = .false.
    logical :: ntrailwaves = .false.
    logical :: nwavesets = .false.
    logical :: packagedata_ifno = .false.
    logical :: cellid = .false.
    logical :: landflag = .false.
    logical :: ivertcon = .false.
    logical :: surfdep = .false.
    logical :: vks = .false.
    logical :: thtr = .false.
    logical :: thts = .false.
    logical :: thti = .false.
    logical :: eps = .false.
    logical :: boundname = .false.
    logical :: ifno = .false.
    logical :: finf = .false.
    logical :: pet = .false.
    logical :: extdp = .false.
    logical :: extwc = .false.
    logical :: ha = .false.
    logical :: hroot = .false.
    logical :: rootact = .false.
    logical :: aux = .false.
  end type GwfUzfParamFoundType

  logical :: gwf_uzf_multi_package = .true.

  character(len=16), parameter :: &
    gwf_uzf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_auxmultname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'AUXMULTNAME', & ! tag name
    'AUXMULTNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of auxiliary variable for multiplier', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
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
    gwfuzf_print_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
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
    gwfuzf_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_wc_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'WC_FILERECORD', & ! tag name
    'WC_FILERECORD', & ! fortran variable
    'RECORD WATER_CONTENT FILEOUT WCFILE', & ! type
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
    gwfuzf_water_content = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'WATER_CONTENT', & ! tag name
    'WATER_CONTENT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'water_content keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_wcfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'WCFILE', & ! tag name
    'WCFILE', & ! fortran variable
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
    gwfuzf_budfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_budget = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_budgetfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_budcsvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_budgetcsv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_budgetcsvfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_pkgconvfilerec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_pkgconv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_pkgconvfname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_mover = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_iconvchk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_NO_FINAL_CHECK', & ! tag name
    'ICONVCHK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not check final convergence', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_simulate_et = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'SIMULATE_ET', & ! tag name
    'SIMULATE_ET', & ! fortran variable
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
    gwfuzf_linear_gwet = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'LINEAR_GWET', & ! tag name
    'LINEAR_GWET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use linear evapotranspiration', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_square_gwet = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'SQUARE_GWET', & ! tag name
    'SQUARE_GWET', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use square evapotranspiration', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_simulate_gwseep = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'SIMULATE_GWSEEP', & ! tag name
    'SIMULATE_GWSEEP', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate seepage', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_unsat_etwc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'UNSAT_ETWC', & ! tag name
    'UNSAT_ETWC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use PET for theta greater than extwc', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_unsat_etae = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'OPTIONS', & ! block
    'UNSAT_ETAE', & ! tag name
    'UNSAT_ETAE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use root potential', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NUZFCELLS', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of UZF cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_ntrailwaves = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NTRAILWAVES', & ! tag name
    'NTRAILWAVES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of trailing waves', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_nwavesets = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NWAVESETS', & ! tag name
    'NWAVESETS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of wave sets', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_packagedata_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IFNO', & ! tag name
    'PACKAGEDATA_IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'uzf id number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_landflag = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'LANDFLAG', & ! tag name
    'LANDFLAG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'land flag', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_ivertcon = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IVERTCON', & ! tag name
    'IVERTCON', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'vertical connection flag', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_surfdep = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'SURFDEP', & ! tag name
    'SURFDEP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'surface depression depth', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_vks = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'VKS', & ! tag name
    'VKS', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'vertical saturated hydraulic conductivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_thtr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'THTR', & ! tag name
    'THTR', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'residual water content', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_thts = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'THTS', & ! tag name
    'THTS', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'saturated water content', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_thti = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'THTI', & ! tag name
    'THTI', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial water content', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_eps = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'EPS', & ! tag name
    'EPS', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'Brooks-Corey exponent', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
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
    gwfuzf_ifno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'IFNO', & ! tag name
    'IFNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'UZF id number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_finf = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'FINF', & ! tag name
    'FINF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'infiltration rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_pet = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'PET', & ! tag name
    'PET', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'potential ET rate', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_extdp = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'EXTDP', & ! tag name
    'EXTDP', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'extinction depth', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_extwc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'EXTWC', & ! tag name
    'EXTWC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'extinction water content', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_ha = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'HA', & ! tag name
    'HA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'air entry potential', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_hroot = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'HROOT', & ! tag name
    'HROOT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'root potential', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_rootact = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'ROOTACT', & ! tag name
    'ROOTACT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'root activity function', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_aux = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
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
    gwf_uzf_param_definitions(*) = &
    [ &
    gwfuzf_auxiliary, &
    gwfuzf_auxmultname, &
    gwfuzf_boundnames, &
    gwfuzf_print_input, &
    gwfuzf_print_flows, &
    gwfuzf_ipakcb, &
    gwfuzf_wc_filerecord, &
    gwfuzf_water_content, &
    gwfuzf_wcfile, &
    gwfuzf_budfilerec, &
    gwfuzf_budget, &
    gwfuzf_fileout, &
    gwfuzf_budgetfile, &
    gwfuzf_budcsvfilerec, &
    gwfuzf_budgetcsv, &
    gwfuzf_budgetcsvfile, &
    gwfuzf_pkgconvfilerec, &
    gwfuzf_pkgconv, &
    gwfuzf_pkgconvfname, &
    gwfuzf_ts_filerecord, &
    gwfuzf_ts6, &
    gwfuzf_filein, &
    gwfuzf_ts6_filename, &
    gwfuzf_obs_filerecord, &
    gwfuzf_obs6, &
    gwfuzf_obs6_filename, &
    gwfuzf_mover, &
    gwfuzf_iconvchk, &
    gwfuzf_simulate_et, &
    gwfuzf_linear_gwet, &
    gwfuzf_square_gwet, &
    gwfuzf_simulate_gwseep, &
    gwfuzf_unsat_etwc, &
    gwfuzf_unsat_etae, &
    gwfuzf_maxbound, &
    gwfuzf_ntrailwaves, &
    gwfuzf_nwavesets, &
    gwfuzf_packagedata_ifno, &
    gwfuzf_cellid, &
    gwfuzf_landflag, &
    gwfuzf_ivertcon, &
    gwfuzf_surfdep, &
    gwfuzf_vks, &
    gwfuzf_thtr, &
    gwfuzf_thts, &
    gwfuzf_thti, &
    gwfuzf_eps, &
    gwfuzf_boundname, &
    gwfuzf_ifno, &
    gwfuzf_finf, &
    gwfuzf_pet, &
    gwfuzf_extdp, &
    gwfuzf_extwc, &
    gwfuzf_ha, &
    gwfuzf_hroot, &
    gwfuzf_rootact, &
    gwfuzf_aux &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfuzf_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IFNO CELLID LANDFLAG IVERTCON SURFDEP VKS THTR THTS '// &
    'THTI EPS BOUNDNAME', & ! type
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
    gwfuzf_perioddata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'UZF', & ! subcomponent
    'PERIOD', & ! block
    'PERIODDATA', & ! tag name
    'PERIODDATA', & ! fortran variable
    'RECARRAY IFNO FINF PET EXTDP EXTWC HA HROOT ROOTACT AUX', & ! type
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
    gwf_uzf_aggregate_definitions(*) = &
    [ &
    gwfuzf_packagedata, &
    gwfuzf_perioddata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_uzf_block_definitions(*) = &
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
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfUzfInputModule

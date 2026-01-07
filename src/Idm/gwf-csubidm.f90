! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfCsubInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_csub_param_definitions
  public gwf_csub_aggregate_definitions
  public gwf_csub_block_definitions
  public GwfCsubParamFoundType
  public gwf_csub_multi_package
  public gwf_csub_subpackages

  type GwfCsubParamFoundType
    logical :: boundnames = .false.
    logical :: print_input = .false.
    logical :: save_flows = .false.
    logical :: gammaw = .false.
    logical :: beta = .false.
    logical :: head_based = .false.
    logical :: precon_head = .false.
    logical :: ndelaycells = .false.
    logical :: icompress = .false.
    logical :: matprop = .false.
    logical :: cell_fraction = .false.
    logical :: interbed_state = .false.
    logical :: precon_stress = .false.
    logical :: delay_head = .false.
    logical :: stress_lag = .false.
    logical :: strainibfr = .false.
    logical :: csvinterbed = .false.
    logical :: fileout = .false.
    logical :: interbedstrainfn = .false.
    logical :: straincgfr = .false.
    logical :: csvcoarse = .false.
    logical :: coarsestrainfn = .false.
    logical :: cmpfr = .false.
    logical :: compaction = .false.
    logical :: cmpfn = .false.
    logical :: cmpelasticfr = .false.
    logical :: cmpelastic = .false.
    logical :: elasticcmpfn = .false.
    logical :: cmpinelasticfr = .false.
    logical :: cmpinelastic = .false.
    logical :: inelasticcmpfn = .false.
    logical :: cmpinterbedfr = .false.
    logical :: cmpinterbed = .false.
    logical :: interbedcmpfn = .false.
    logical :: cmpcoarsefr = .false.
    logical :: cmpcoarse = .false.
    logical :: cmpcoarsefn = .false.
    logical :: zdispfr = .false.
    logical :: zdisplacement = .false.
    logical :: zdispfn = .false.
    logical :: pkgconvergefr = .false.
    logical :: pkgconverge = .false.
    logical :: pkgconvergefn = .false.
    logical :: ts_filerecord = .false.
    logical :: ts6 = .false.
    logical :: filein = .false.
    logical :: ts6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: ninterbeds = .false.
    logical :: maxbound = .false.
    logical :: cg_ske_cr = .false.
    logical :: cg_theta = .false.
    logical :: sgm = .false.
    logical :: sgs = .false.
    logical :: icsubno = .false.
    logical :: cellid_pkgdata = .false.
    logical :: cdelay = .false.
    logical :: pcs0 = .false.
    logical :: thick_frac = .false.
    logical :: rnb = .false.
    logical :: ssv_cc = .false.
    logical :: sse_cr = .false.
    logical :: theta = .false.
    logical :: kv = .false.
    logical :: h0 = .false.
    logical :: boundname = .false.
    logical :: cellid = .false.
    logical :: sig0 = .false.
  end type GwfCsubParamFoundType

  logical :: gwf_csub_multi_package = .false.

  character(len=16), parameter :: &
    gwf_csub_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_save_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save CSUB flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_gammaw = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'GAMMAW', & ! tag name
    'GAMMAW', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'unit weight of water', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_beta = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'BETA', & ! tag name
    'BETA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'compressibility of water', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_head_based = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'HEAD_BASED', & ! tag name
    'HEAD_BASED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate the head-based formulation will be used', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_precon_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'INITIAL_PRECONSOLIDATION_HEAD', & ! tag name
    'PRECON_HEAD', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that preconsolidation heads will be specified', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_ndelaycells = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'NDELAYCELLS', & ! tag name
    'NDELAYCELLS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of interbed cell nodes', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_icompress = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPRESSION_INDICES', & ! tag name
    'ICOMPRESS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate CR and CC are read instead of SSE and SSV', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_matprop = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'UPDATE_MATERIAL_PROPERTIES', & ! tag name
    'MATPROP', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate material properties can change during the&
     & simulations', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cell_fraction = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'CELL_FRACTION', & ! tag name
    'CELL_FRACTION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate cell fraction interbed thickness', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_interbed_state = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'SPECIFIED_INITIAL_INTERBED_STATE', & ! tag name
    'INTERBED_STATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that absolute initial states will be specified', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_precon_stress = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'SPECIFIED_INITIAL_PRECONSOLIDATION_STRESS', & ! tag name
    'PRECON_STRESS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that absolute initial preconsolidation stresses&
     & (head) will be specified', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_delay_head = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'SPECIFIED_INITIAL_DELAY_HEAD', & ! tag name
    'DELAY_HEAD', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that absolute initial delay bed heads will be&
     & specified', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_stress_lag = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'EFFECTIVE_STRESS_LAG', & ! tag name
    'STRESS_LAG', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate that specific storage will be calculate using the&
     & effective stress from the previous time step', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_strainibfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'STRAINIB_FILERECORD', & ! tag name
    'STRAINIBFR', & ! fortran variable
    'RECORD STRAIN_CSV_INTERBED FILEOUT INTERBEDSTRAIN_FILENAME', & ! type
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
    gwfcsub_csvinterbed = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'STRAIN_CSV_INTERBED', & ! tag name
    'CSVINTERBED', & ! fortran variable
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
    gwfcsub_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_interbedstrainfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'INTERBEDSTRAIN_FILENAME', & ! tag name
    'INTERBEDSTRAINFN', & ! fortran variable
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
    gwfcsub_straincgfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'STRAINCG_FILERECORD', & ! tag name
    'STRAINCGFR', & ! fortran variable
    'RECORD STRAIN_CSV_COARSE FILEOUT COARSESTRAIN_FILENAME', & ! type
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
    gwfcsub_csvcoarse = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'STRAIN_CSV_COARSE', & ! tag name
    'CSVCOARSE', & ! fortran variable
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
    gwfcsub_coarsestrainfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COARSESTRAIN_FILENAME', & ! tag name
    'COARSESTRAINFN', & ! fortran variable
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
    gwfcsub_cmpfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_FILERECORD', & ! tag name
    'CMPFR', & ! fortran variable
    'RECORD COMPACTION FILEOUT COMPACTION_FILENAME', & ! type
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
    gwfcsub_compaction = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION', & ! tag name
    'COMPACTION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'compaction keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cmpfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_FILENAME', & ! tag name
    'CMPFN', & ! fortran variable
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
    gwfcsub_cmpelasticfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_ELASTIC_FILERECORD', & ! tag name
    'CMPELASTICFR', & ! fortran variable
    'RECORD COMPACTION_ELASTIC FILEOUT ELASTIC_COMPACTION_FILENAME', & ! type
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
    gwfcsub_cmpelastic = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_ELASTIC', & ! tag name
    'CMPELASTIC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'elastic interbed compaction keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_elasticcmpfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'ELASTIC_COMPACTION_FILENAME', & ! tag name
    'ELASTICCMPFN', & ! fortran variable
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
    gwfcsub_cmpinelasticfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_INELASTIC_FILERECORD', & ! tag name
    'CMPINELASTICFR', & ! fortran variable
    'RECORD COMPACTION_INELASTIC FILEOUT INELASTIC_COMPACTION_FILENAME', & ! type
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
    gwfcsub_cmpinelastic = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_INELASTIC', & ! tag name
    'CMPINELASTIC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'inelastic interbed compaction keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_inelasticcmpfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'INELASTIC_COMPACTION_FILENAME', & ! tag name
    'INELASTICCMPFN', & ! fortran variable
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
    gwfcsub_cmpinterbedfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_INTERBED_FILERECORD', & ! tag name
    'CMPINTERBEDFR', & ! fortran variable
    'RECORD COMPACTION_INTERBED FILEOUT INTERBED_COMPACTION_FILENAME', & ! type
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
    gwfcsub_cmpinterbed = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_INTERBED', & ! tag name
    'CMPINTERBED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'interbed compaction keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_interbedcmpfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'INTERBED_COMPACTION_FILENAME', & ! tag name
    'INTERBEDCMPFN', & ! fortran variable
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
    gwfcsub_cmpcoarsefr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_COARSE_FILERECORD', & ! tag name
    'CMPCOARSEFR', & ! fortran variable
    'RECORD COMPACTION_COARSE FILEOUT COARSE_COMPACTION_FILENAME', & ! type
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
    gwfcsub_cmpcoarse = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COMPACTION_COARSE', & ! tag name
    'CMPCOARSE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'coarse compaction keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cmpcoarsefn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'COARSE_COMPACTION_FILENAME', & ! tag name
    'CMPCOARSEFN', & ! fortran variable
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
    gwfcsub_zdispfr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'ZDISPLACEMENT_FILERECORD', & ! tag name
    'ZDISPFR', & ! fortran variable
    'RECORD ZDISPLACEMENT FILEOUT ZDISPLACEMENT_FILENAME', & ! type
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
    gwfcsub_zdisplacement = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'ZDISPLACEMENT', & ! tag name
    'ZDISPLACEMENT', & ! fortran variable
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
    gwfcsub_zdispfn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'ZDISPLACEMENT_FILENAME', & ! tag name
    'ZDISPFN', & ! fortran variable
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
    gwfcsub_pkgconvergefr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILERECORD', & ! tag name
    'PKGCONVERGEFR', & ! fortran variable
    'RECORD PACKAGE_CONVERGENCE FILEOUT PACKAGE_CONVERGENCE_FILENAME', & ! type
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
    gwfcsub_pkgconverge = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE', & ! tag name
    'PKGCONVERGE', & ! fortran variable
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
    gwfcsub_pkgconvergefn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'OPTIONS', & ! block
    'PACKAGE_CONVERGENCE_FILENAME', & ! tag name
    'PKGCONVERGEFN', & ! fortran variable
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
    gwfcsub_ts_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_ts6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_ts6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_obs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_obs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_obs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_ninterbeds = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'DIMENSIONS', & ! block
    'NINTERBEDS', & ! tag name
    'NINTERBEDS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of CSUB interbed systems', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXSIG0', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of stress offset cells', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cg_ske_cr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'GRIDDATA', & ! block
    'CG_SKE_CR', & ! tag name
    'CG_SKE_CR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'elastic coarse specific storage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cg_theta = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'GRIDDATA', & ! block
    'CG_THETA', & ! tag name
    'CG_THETA', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'initial coarse-grained material porosity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_sgm = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'GRIDDATA', & ! block
    'SGM', & ! tag name
    'SGM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'specific gravity of moist sediments', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_sgs = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'GRIDDATA', & ! block
    'SGS', & ! tag name
    'SGS', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'specific gravity of saturated sediments', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_icsubno = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'ICSUBNO', & ! tag name
    'ICSUBNO', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'CSUB id number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_cellid_pkgdata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'CELLID', & ! tag name
    'CELLID_PKGDATA', & ! fortran variable
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
    gwfcsub_cdelay = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'CDELAY', & ! tag name
    'CDELAY', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'delay type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_pcs0 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PCS0', & ! tag name
    'PCS0', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial stress', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_thick_frac = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'THICK_FRAC', & ! tag name
    'THICK_FRAC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'interbed thickness or cell fraction', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_rnb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'RNB', & ! tag name
    'RNB', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'delay interbed material factor', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_ssv_cc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'SSV_CC', & ! tag name
    'SSV_CC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial interbed inelastic specific storage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_sse_cr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'SSE_CR', & ! tag name
    'SSE_CR', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial interbed elastic specific storage', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_theta = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'THETA', & ! tag name
    'THETA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial interbed porosity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_kv = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'KV', & ! tag name
    'KV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'delay interbed vertical hydraulic conductivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_h0 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'H0', & ! tag name
    'H0', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'initial delay interbed head', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
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
    gwfcsub_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PERIOD', & ! block
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
    gwfcsub_sig0 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PERIOD', & ! block
    'SIG0', & ! tag name
    'SIG0', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'well stress offset', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .true. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_csub_param_definitions(*) = &
    [ &
    gwfcsub_boundnames, &
    gwfcsub_print_input, &
    gwfcsub_save_flows, &
    gwfcsub_gammaw, &
    gwfcsub_beta, &
    gwfcsub_head_based, &
    gwfcsub_precon_head, &
    gwfcsub_ndelaycells, &
    gwfcsub_icompress, &
    gwfcsub_matprop, &
    gwfcsub_cell_fraction, &
    gwfcsub_interbed_state, &
    gwfcsub_precon_stress, &
    gwfcsub_delay_head, &
    gwfcsub_stress_lag, &
    gwfcsub_strainibfr, &
    gwfcsub_csvinterbed, &
    gwfcsub_fileout, &
    gwfcsub_interbedstrainfn, &
    gwfcsub_straincgfr, &
    gwfcsub_csvcoarse, &
    gwfcsub_coarsestrainfn, &
    gwfcsub_cmpfr, &
    gwfcsub_compaction, &
    gwfcsub_cmpfn, &
    gwfcsub_cmpelasticfr, &
    gwfcsub_cmpelastic, &
    gwfcsub_elasticcmpfn, &
    gwfcsub_cmpinelasticfr, &
    gwfcsub_cmpinelastic, &
    gwfcsub_inelasticcmpfn, &
    gwfcsub_cmpinterbedfr, &
    gwfcsub_cmpinterbed, &
    gwfcsub_interbedcmpfn, &
    gwfcsub_cmpcoarsefr, &
    gwfcsub_cmpcoarse, &
    gwfcsub_cmpcoarsefn, &
    gwfcsub_zdispfr, &
    gwfcsub_zdisplacement, &
    gwfcsub_zdispfn, &
    gwfcsub_pkgconvergefr, &
    gwfcsub_pkgconverge, &
    gwfcsub_pkgconvergefn, &
    gwfcsub_ts_filerecord, &
    gwfcsub_ts6, &
    gwfcsub_filein, &
    gwfcsub_ts6_filename, &
    gwfcsub_obs_filerecord, &
    gwfcsub_obs6, &
    gwfcsub_obs6_filename, &
    gwfcsub_ninterbeds, &
    gwfcsub_maxbound, &
    gwfcsub_cg_ske_cr, &
    gwfcsub_cg_theta, &
    gwfcsub_sgm, &
    gwfcsub_sgs, &
    gwfcsub_icsubno, &
    gwfcsub_cellid_pkgdata, &
    gwfcsub_cdelay, &
    gwfcsub_pcs0, &
    gwfcsub_thick_frac, &
    gwfcsub_rnb, &
    gwfcsub_ssv_cc, &
    gwfcsub_sse_cr, &
    gwfcsub_theta, &
    gwfcsub_kv, &
    gwfcsub_h0, &
    gwfcsub_boundname, &
    gwfcsub_cellid, &
    gwfcsub_sig0 &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY ICSUBNO CELLID CDELAY PCS0 THICK_FRAC RNB SSV_CC SSE_CR THETA KV H0 BOUNDNAME', & ! type
    'NINTERBEDS', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfcsub_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'CSUB', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID SIG0', & ! type
    'MAXSIG0', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_csub_aggregate_definitions(*) = &
    [ &
    gwfcsub_packagedata, &
    gwfcsub_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_csub_block_definitions(*) = &
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
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGEDATA', & ! blockname
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

end module GwfCsubInputModule

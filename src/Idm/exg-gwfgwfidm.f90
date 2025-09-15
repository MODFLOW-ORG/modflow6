! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgGwfgwfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_gwfgwf_param_definitions
  public exg_gwfgwf_aggregate_definitions
  public exg_gwfgwf_block_definitions
  public ExgGwfgwfParamFoundType
  public exg_gwfgwf_multi_package
  public exg_gwfgwf_subpackages

  type ExgGwfgwfParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: cell_averaging = .false.
    logical :: cvoptions = .false.
    logical :: variablecv = .false.
    logical :: dewatered = .false.
    logical :: newton = .false.
    logical :: xt3d = .false.
    logical :: gnc_filerecord = .false.
    logical :: filein = .false.
    logical :: gnc6 = .false.
    logical :: gnc6_filename = .false.
    logical :: mvr_filerecord = .false.
    logical :: mvr6 = .false.
    logical :: mvr6_filename = .false.
    logical :: obs_filerecord = .false.
    logical :: obs6 = .false.
    logical :: obs6_filename = .false.
    logical :: dev_ifmod_on = .false.
    logical :: nexg = .false.
    logical :: cellidm1 = .false.
    logical :: cellidm2 = .false.
    logical :: ihc = .false.
    logical :: cl1 = .false.
    logical :: cl2 = .false.
    logical :: hwva = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type ExgGwfgwfParamFoundType

  logical :: exg_gwfgwf_multi_package = .true.

  character(len=16), parameter :: &
    exg_gwfgwf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_auxiliary = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_boundnames = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_iprpak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print input to list file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_iprflow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print gwfgwf flows to list file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_ipakcb = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save GWFGWF flows', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cell_averaging = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'CELL_AVERAGING', & ! tag name
    'CELL_AVERAGING', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'conductance weighting option', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cvoptions = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'CVOPTIONS', & ! tag name
    'CVOPTIONS', & ! fortran variable
    'RECORD VARIABLECV DEWATERED', & ! type
    '', & ! shape
    'vertical conductance options', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_variablecv = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'VARIABLECV', & ! tag name
    'VARIABLECV', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate VARIABLECV option', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_dewatered = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'DEWATERED', & ! tag name
    'DEWATERED', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate DEWATERED option', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_newton = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTON', & ! tag name
    'NEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_xt3d = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'XT3D', & ! tag name
    'XT3D', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate XT3D', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_gnc_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'GNC_FILERECORD', & ! tag name
    'GNC_FILERECORD', & ! fortran variable
    'RECORD GNC6 FILEIN GNC6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_filein = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_gnc6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'GNC6', & ! tag name
    'GNC6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'gnc6 keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_gnc6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'GNC6_FILENAME', & ! tag name
    'GNC6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'gnc6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_mvr_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'MVR_FILERECORD', & ! tag name
    'MVR_FILERECORD', & ! fortran variable
    'RECORD MVR6 FILEIN MVR6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_mvr6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'MVR6', & ! tag name
    'MVR6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_mvr6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'MVR6_FILENAME', & ! tag name
    'MVR6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'mvr6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_obs_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS_FILERECORD', & ! tag name
    'OBS_FILERECORD', & ! fortran variable
    'RECORD OBS6 FILEIN OBS6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_obs6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6', & ! tag name
    'OBS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'obs keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_obs6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'OBS6_FILENAME', & ! tag name
    'OBS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'obs6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_dev_ifmod_on = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_INTERFACEMODEL_ON', & ! tag name
    'DEV_IFMOD_ON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate interface model on exchange', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'DIMENSIONS', & ! block
    'NEXG', & ! tag name
    'NEXG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of exchanges', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM1', & ! tag name
    'CELLIDM1', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of first cell', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM2', & ! tag name
    'CELLIDM2', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of second cell', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_ihc = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'IHC', & ! tag name
    'IHC', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'integer flag for connection type', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cl1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CL1', & ! tag name
    'CL1', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection distance', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_cl2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CL2', & ! tag name
    'CL2', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection distance', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_hwva = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'horizontal cell width or area for vertical flow', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_auxvar = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_boundname = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'exchange boundname', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_gwfgwf_param_definitions(*) = &
    [ &
    exggwfgwf_auxiliary, &
    exggwfgwf_boundnames, &
    exggwfgwf_iprpak, &
    exggwfgwf_iprflow, &
    exggwfgwf_ipakcb, &
    exggwfgwf_cell_averaging, &
    exggwfgwf_cvoptions, &
    exggwfgwf_variablecv, &
    exggwfgwf_dewatered, &
    exggwfgwf_newton, &
    exggwfgwf_xt3d, &
    exggwfgwf_gnc_filerecord, &
    exggwfgwf_filein, &
    exggwfgwf_gnc6, &
    exggwfgwf_gnc6_filename, &
    exggwfgwf_mvr_filerecord, &
    exggwfgwf_mvr6, &
    exggwfgwf_mvr6_filename, &
    exggwfgwf_obs_filerecord, &
    exggwfgwf_obs6, &
    exggwfgwf_obs6_filename, &
    exggwfgwf_dev_ifmod_on, &
    exggwfgwf_nexg, &
    exggwfgwf_cellidm1, &
    exggwfgwf_cellidm2, &
    exggwfgwf_ihc, &
    exggwfgwf_cl1, &
    exggwfgwf_cl2, &
    exggwfgwf_hwva, &
    exggwfgwf_auxvar, &
    exggwfgwf_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwfgwf_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWFGWF', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'EXCHANGEDATA', & ! tag name
    'EXCHANGEDATA', & ! fortran variable
    'RECARRAY CELLIDM1 CELLIDM2 IHC CL1 CL2 HWVA AUX BOUNDNAME', & ! type
    'NEXG', & ! shape
    'exchange data', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_gwfgwf_aggregate_definitions(*) = &
    [ &
    exggwfgwf_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_gwfgwf_block_definitions(*) = &
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
    'EXCHANGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ExgGwfgwfInputModule

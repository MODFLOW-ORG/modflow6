! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgGwegweInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_gwegwe_param_definitions
  public exg_gwegwe_aggregate_definitions
  public exg_gwegwe_block_definitions
  public ExgGwegweParamFoundType
  public exg_gwegwe_multi_package
  public exg_gwegwe_subpackages

  type ExgGwegweParamFoundType
    logical :: gwfmodelname1 = .false.
    logical :: gwfmodelname2 = .false.
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: adv_scheme = .false.
    logical :: cnd_xt3d_off = .false.
    logical :: cnd_xt3d_rhs = .false.
    logical :: filein = .false.
    logical :: mve_filerecord = .false.
    logical :: mve6 = .false.
    logical :: mve6_filename = .false.
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
  end type ExgGwegweParamFoundType

  logical :: exg_gwegwe_multi_package = .true.

  character(len=16), parameter :: &
    exg_gwegwe_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_gwfmodelname1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME1', & ! tag name
    'GWFMODELNAME1', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of first corresponding GWF Model', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_gwfmodelname2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME2', & ! tag name
    'GWFMODELNAME2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of second corresponding GWF Model', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_auxiliary = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_boundnames = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_iprpak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_iprflow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_ipakcb = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_adv_scheme = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'ADV_SCHEME', & ! tag name
    'ADV_SCHEME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'advective scheme', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_cnd_xt3d_off = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'CND_XT3D_OFF', & ! tag name
    'CND_XT3D_OFF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'deactivate xt3d', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_cnd_xt3d_rhs = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'CND_XT3D_RHS', & ! tag name
    'CND_XT3D_RHS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'xt3d on right-hand side', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_filein = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_mve_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'MVE_FILERECORD', & ! tag name
    'MVE_FILERECORD', & ! fortran variable
    'RECORD MVE6 FILEIN MVE6_FILENAME', & ! type
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
    exggwegwe_mve6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'MVE6', & ! tag name
    'MVE6', & ! fortran variable
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
    exggwegwe_mve6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
    'OPTIONS', & ! block
    'MVE6_FILENAME', & ! tag name
    'MVE6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'mve6 input filename', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_obs_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_obs6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_obs6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_dev_ifmod_on = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_ihc = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_cl1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_cl2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_hwva = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_auxvar = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exggwegwe_boundname = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exg_gwegwe_param_definitions(*) = &
    [ &
    exggwegwe_gwfmodelname1, &
    exggwegwe_gwfmodelname2, &
    exggwegwe_auxiliary, &
    exggwegwe_boundnames, &
    exggwegwe_iprpak, &
    exggwegwe_iprflow, &
    exggwegwe_ipakcb, &
    exggwegwe_adv_scheme, &
    exggwegwe_cnd_xt3d_off, &
    exggwegwe_cnd_xt3d_rhs, &
    exggwegwe_filein, &
    exggwegwe_mve_filerecord, &
    exggwegwe_mve6, &
    exggwegwe_mve6_filename, &
    exggwegwe_obs_filerecord, &
    exggwegwe_obs6, &
    exggwegwe_obs6_filename, &
    exggwegwe_dev_ifmod_on, &
    exggwegwe_nexg, &
    exggwegwe_cellidm1, &
    exggwegwe_cellidm2, &
    exggwegwe_ihc, &
    exggwegwe_cl1, &
    exggwegwe_cl2, &
    exggwegwe_hwva, &
    exggwegwe_auxvar, &
    exggwegwe_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwegwe_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWEGWE', & ! subcomponent
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
    exg_gwegwe_aggregate_definitions(*) = &
    [ &
    exggwegwe_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_gwegwe_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .true., & ! required
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

end module ExgGwegweInputModule

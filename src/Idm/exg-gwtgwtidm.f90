! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgGwtgwtInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_gwtgwt_param_definitions
  public exg_gwtgwt_aggregate_definitions
  public exg_gwtgwt_block_definitions
  public ExgGwtgwtParamFoundType
  public exg_gwtgwt_multi_package
  public exg_gwtgwt_subpackages

  type ExgGwtgwtParamFoundType
    logical :: gwfmodelname1 = .false.
    logical :: gwfmodelname2 = .false.
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: adv_scheme = .false.
    logical :: dsp_xt3d_off = .false.
    logical :: dsp_xt3d_rhs = .false.
    logical :: filein = .false.
    logical :: mvt_filerecord = .false.
    logical :: mvt6 = .false.
    logical :: mvt6_filename = .false.
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
  end type ExgGwtgwtParamFoundType

  logical :: exg_gwtgwt_multi_package = .true.

  character(len=16), parameter :: &
    exg_gwtgwt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_gwfmodelname1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME1', & ! tag name
    'GWFMODELNAME1', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of first corresponding GWF Model', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_gwfmodelname2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME2', & ! tag name
    'GWFMODELNAME2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'keyword to specify name of second corresponding GWF Model', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_auxiliary = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_boundnames = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_iprpak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print input to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_iprflow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print gwfgwf flows to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_ipakcb = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save GWFGWF flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_adv_scheme = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'ADV_SCHEME', & ! tag name
    'ADV_SCHEME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'advective scheme', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_dsp_xt3d_off = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'DSP_XT3D_OFF', & ! tag name
    'DSP_XT3D_OFF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'deactivate xt3d', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_dsp_xt3d_rhs = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'DSP_XT3D_RHS', & ! tag name
    'DSP_XT3D_RHS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'xt3d on right-hand side', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_filein = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_mvt_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'MVT_FILERECORD', & ! tag name
    'MVT_FILERECORD', & ! fortran variable
    'RECORD MVT6 FILEIN MVT6_FILENAME', & ! type
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
    exggwtgwt_mvt6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'MVT6', & ! tag name
    'MVT6', & ! fortran variable
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
    exggwtgwt_mvt6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'MVT6_FILENAME', & ! tag name
    'MVT6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'mvt6 input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_obs_filerecord = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_obs6 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_obs6_filename = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
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
    exggwtgwt_dev_ifmod_on = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_INTERFACEMODEL_ON', & ! tag name
    'DEV_IFMOD_ON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'activate interface model on exchange', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'DIMENSIONS', & ! block
    'NEXG', & ! tag name
    'NEXG', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of exchanges', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM1', & ! tag name
    'CELLIDM1', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of first cell', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CELLIDM2', & ! tag name
    'CELLIDM2', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cellid of second cell', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_ihc = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'IHC', & ! tag name
    'IHC', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'integer flag for connection type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_cl1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CL1', & ! tag name
    'CL1', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection distance', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_cl2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'CL2', & ! tag name
    'CL2', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'connection distance', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_hwva = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'horizontal cell width or area for vertical flow', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_auxvar = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_boundname = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'exchange boundname', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_gwtgwt_param_definitions(*) = &
    [ &
    exggwtgwt_gwfmodelname1, &
    exggwtgwt_gwfmodelname2, &
    exggwtgwt_auxiliary, &
    exggwtgwt_boundnames, &
    exggwtgwt_iprpak, &
    exggwtgwt_iprflow, &
    exggwtgwt_ipakcb, &
    exggwtgwt_adv_scheme, &
    exggwtgwt_dsp_xt3d_off, &
    exggwtgwt_dsp_xt3d_rhs, &
    exggwtgwt_filein, &
    exggwtgwt_mvt_filerecord, &
    exggwtgwt_mvt6, &
    exggwtgwt_mvt6_filename, &
    exggwtgwt_obs_filerecord, &
    exggwtgwt_obs6, &
    exggwtgwt_obs6_filename, &
    exggwtgwt_dev_ifmod_on, &
    exggwtgwt_nexg, &
    exggwtgwt_cellidm1, &
    exggwtgwt_cellidm2, &
    exggwtgwt_ihc, &
    exggwtgwt_cl1, &
    exggwtgwt_cl2, &
    exggwtgwt_hwva, &
    exggwtgwt_auxvar, &
    exggwtgwt_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    exggwtgwt_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'GWTGWT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'EXCHANGEDATA', & ! tag name
    'EXCHANGEDATA', & ! fortran variable
    'RECARRAY CELLIDM1 CELLIDM2 IHC CL1 CL2 HWVA AUX BOUNDNAME', & ! type
    'NEXG', & ! shape
    'exchange data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exg_gwtgwt_aggregate_definitions(*) = &
    [ &
    exggwtgwt_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_gwtgwt_block_definitions(*) = &
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

end module ExgGwtgwtInputModule

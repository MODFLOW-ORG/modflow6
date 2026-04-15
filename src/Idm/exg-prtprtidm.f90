! ** Do Not Modify! MODFLOW 6 system generated file. **
module ExgPrtprtInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public exg_prtprt_param_definitions
  public exg_prtprt_aggregate_definitions
  public exg_prtprt_block_definitions
  public ExgPrtprtParamFoundType
  public exg_prtprt_multi_package
  public exg_prtprt_subpackages

  type ExgPrtprtParamFoundType
    logical :: gwfmodelname1 = .false.
    logical :: gwfmodelname2 = .false.
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: nexg = .false.
    logical :: cellidm1 = .false.
    logical :: cellidm2 = .false.
    logical :: ihc = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type ExgPrtprtParamFoundType

  logical :: exg_prtprt_multi_package = .true.

  character(len=16), parameter :: &
    exg_prtprt_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_gwfmodelname1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME1', & ! tag name
    'GWFMODELNAME1', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of gwf model for prt model 1', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_gwfmodelname2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
    'OPTIONS', & ! block
    'GWFMODELNAME2', & ! tag name
    'GWFMODELNAME2', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of gwf model for prt model 2', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_auxiliary = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_boundnames = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_iprpak = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_iprflow = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to print prtprt flows to list file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_ipakcb = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save PRTPRT flows', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_nexg = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_cellidm1 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_cellidm2 = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_ihc = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_auxvar = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exgprtprt_boundname = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
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
    exg_prtprt_param_definitions(*) = &
    [ &
    exgprtprt_gwfmodelname1, &
    exgprtprt_gwfmodelname2, &
    exgprtprt_auxiliary, &
    exgprtprt_boundnames, &
    exgprtprt_iprpak, &
    exgprtprt_iprflow, &
    exgprtprt_ipakcb, &
    exgprtprt_nexg, &
    exgprtprt_cellidm1, &
    exgprtprt_cellidm2, &
    exgprtprt_ihc, &
    exgprtprt_auxvar, &
    exgprtprt_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    exgprtprt_exchangedata = InputParamDefinitionType &
    ( &
    'EXG', & ! component
    'PRTPRT', & ! subcomponent
    'EXCHANGEDATA', & ! block
    'EXCHANGEDATA', & ! tag name
    'EXCHANGEDATA', & ! fortran variable
    'RECARRAY CELLIDM1 CELLIDM2 IHC AUX BOUNDNAME', & ! type
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
    exg_prtprt_aggregate_definitions(*) = &
    [ &
    exgprtprt_exchangedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    exg_prtprt_block_definitions(*) = &
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

end module ExgPrtprtInputModule

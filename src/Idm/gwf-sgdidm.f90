! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfSgdInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_sgd_param_definitions
  public gwf_sgd_aggregate_definitions
  public gwf_sgd_block_definitions
  public GwfSgdParamFoundType
  public gwf_sgd_multi_package
  public gwf_sgd_is_advanced
  public gwf_sgd_subpackages

  type GwfSgdParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: gradx = .false.
    logical :: grady = .false.
    logical :: gradz = .false.
    logical :: hwva = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwfSgdParamFoundType

  logical :: gwf_sgd_multi_package = .true.
  logical :: gwf_sgd_is_advanced = .false.

  character(len=16), parameter :: &
    gwf_sgd_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
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
    gwfsgd_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
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
    gwfsgd_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
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
    gwfsgd_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print boundary flow to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save boundary flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of specified gradient cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
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
    gwfsgd_gradx = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'GRADX', & ! tag name
    'GRADX', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'specified gradient x component', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_grady = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'GRADY', & ! tag name
    'GRADY', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'specified gradient y component', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_gradz = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'GRADZ', & ! tag name
    'GRADZ', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'specified gradient z component', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_hwva = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'boundary face area', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
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
    gwfsgd_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'specified gradient boundary name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_sgd_param_definitions(*) = &
    [ &
    gwfsgd_auxiliary, &
    gwfsgd_boundnames, &
    gwfsgd_iprpak, &
    gwfsgd_iprflow, &
    gwfsgd_ipakcb, &
    gwfsgd_maxbound, &
    gwfsgd_cellid, &
    gwfsgd_gradx, &
    gwfsgd_grady, &
    gwfsgd_gradz, &
    gwfsgd_hwva, &
    gwfsgd_auxvar, &
    gwfsgd_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfsgd_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SGD', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID GRADX GRADY GRADZ HWVA AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_sgd_aggregate_definitions(*) = &
    [ &
    gwfsgd_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_sgd_block_definitions(*) = &
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
    'PERIOD', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfSgdInputModule

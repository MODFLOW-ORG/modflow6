! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfSpgInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_spg_param_definitions
  public gwf_spg_aggregate_definitions
  public gwf_spg_block_definitions
  public GwfSpgParamFoundType
  public gwf_spg_multi_package
  public gwf_spg_is_advanced
  public gwf_spg_subpackages

  type GwfSpgParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: ihc = .false.
    logical :: cl1 = .false.
    logical :: hwva = .false.
    logical :: angldegx = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwfSpgParamFoundType

  logical :: gwf_spg_multi_package = .true.
  logical :: gwf_spg_is_advanced = .false.

  character(len=16), parameter :: &
    gwf_spg_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfspg_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
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
    gwfspg_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
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
    gwfspg_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
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
    gwfspg_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print seepage rates to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save seepage flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of seepage cells', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
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
    gwfspg_ihc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'IHC', & ! tag name
    'IHC', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'horizontal connection flag', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_cl1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'CL1', & ! tag name
    'CL1', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'distance to seepage face', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_hwva = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'seepage face area', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_angldegx = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'ANGLDEGX', & ! tag name
    'ANGLDEGX', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'face normal angle with x axis', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspg_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
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
    gwfspg_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'seepage boundary name', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_spg_param_definitions(*) = &
    [ &
    gwfspg_auxiliary, &
    gwfspg_boundnames, &
    gwfspg_iprpak, &
    gwfspg_iprflow, &
    gwfspg_ipakcb, &
    gwfspg_maxbound, &
    gwfspg_cellid, &
    gwfspg_ihc, &
    gwfspg_cl1, &
    gwfspg_hwva, &
    gwfspg_angldegx, &
    gwfspg_auxvar, &
    gwfspg_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfspg_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPG', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID IHC CL1 HWVA ANGLDEGX AUX BOUNDNAME', & ! type
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
    gwf_spg_aggregate_definitions(*) = &
    [ &
    gwfspg_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_spg_block_definitions(*) = &
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

end module GwfSpgInputModule

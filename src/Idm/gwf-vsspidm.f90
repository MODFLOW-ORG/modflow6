! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfVsspInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_vssp_param_definitions
  public gwf_vssp_aggregate_definitions
  public gwf_vssp_block_definitions
  public GwfVsspParamFoundType
  public gwf_vssp_multi_package
  public gwf_vssp_is_advanced
  public gwf_vssp_subpackages

  type GwfVsspParamFoundType
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
  end type GwfVsspParamFoundType

  logical :: gwf_vssp_multi_package = .true.
  logical :: gwf_vssp_is_advanced = .false.

  character(len=16), parameter :: &
    gwf_vssp_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvssp_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_ihc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_cl1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_hwva = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_angldegx = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwfvssp_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwf_vssp_param_definitions(*) = &
    [ &
    gwfvssp_auxiliary, &
    gwfvssp_boundnames, &
    gwfvssp_iprpak, &
    gwfvssp_iprflow, &
    gwfvssp_ipakcb, &
    gwfvssp_maxbound, &
    gwfvssp_cellid, &
    gwfvssp_ihc, &
    gwfvssp_cl1, &
    gwfvssp_hwva, &
    gwfvssp_angldegx, &
    gwfvssp_auxvar, &
    gwfvssp_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvssp_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSSP', & ! subcomponent
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
    gwf_vssp_aggregate_definitions(*) = &
    [ &
    gwfvssp_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_vssp_block_definitions(*) = &
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

end module GwfVsspInputModule

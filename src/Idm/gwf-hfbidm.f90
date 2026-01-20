! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfHfbInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_hfb_param_definitions
  public gwf_hfb_aggregate_definitions
  public gwf_hfb_block_definitions
  public GwfHfbParamFoundType
  public gwf_hfb_multi_package
  public gwf_hfb_subpackages

  type GwfHfbParamFoundType
    logical :: print_input = .false.
    logical :: maxbound = .false.
    logical :: cellid1 = .false.
    logical :: cellid2 = .false.
    logical :: hydchr = .false.
  end type GwfHfbParamFoundType

  logical :: gwf_hfb_multi_package = .false.

  character(len=16), parameter :: &
    gwf_hfb_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'model print input to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXHFB', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of barriers', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_cellid1 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'PERIOD', & ! block
    'CELLID1', & ! tag name
    'CELLID1', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'first cell adjacent to barrier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_cellid2 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'PERIOD', & ! block
    'CELLID2', & ! tag name
    'CELLID2', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'second cell adjacent to barrier', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_hydchr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'PERIOD', & ! block
    'HYDCHR', & ! tag name
    'HYDCHR', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'barrier hydraulic characteristic', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_hfb_param_definitions(*) = &
    [ &
    gwfhfb_print_input, &
    gwfhfb_maxbound, &
    gwfhfb_cellid1, &
    gwfhfb_cellid2, &
    gwfhfb_hydchr &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfhfb_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'HFB', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID1 CELLID2 HYDCHR', & ! type
    'MAXHFB', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_hfb_aggregate_definitions(*) = &
    [ &
    gwfhfb_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_hfb_block_definitions(*) = &
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

end module GwfHfbInputModule

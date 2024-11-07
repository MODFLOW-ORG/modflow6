! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfSpfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_spf_param_definitions
  public gwf_spf_aggregate_definitions
  public gwf_spf_block_definitions
  public GwfSpfParamFoundType
  public gwf_spf_multi_package
  public gwf_spf_subpackages

  type GwfSpfParamFoundType
    logical :: auxiliary = .false.
    logical :: boundnames = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: ipakcb = .false.
    logical :: some_option = .false.
    logical :: maxbound = .false.
    logical :: cellid = .false.
    logical :: dist = .false.
    logical :: area = .false.
    logical :: auxvar = .false.
    logical :: boundname = .false.
  end type GwfSpfParamFoundType

  logical :: gwf_spf_multi_package = .true.

  character(len=16), parameter :: &
    gwf_spf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfspf_auxiliary = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'AUXILIARY', & ! tag name
    'AUXILIARY', & ! fortran variable
    'STRING', & ! type
    'NAUX', & ! shape
    'keyword to specify aux variables', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_boundnames = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'BOUNDNAMES', & ! tag name
    'BOUNDNAMES', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_iprpak = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_iprflow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'IPRFLOW', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print seepage rates to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save seepage flows to budget file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_some_option = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'OPTIONS', & ! block
    'SOME_OPTION', & ! tag name
    'SOME_OPTION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'some option', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_maxbound = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'DIMENSIONS', & ! block
    'MAXBOUND', & ! tag name
    'MAXBOUND', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'maximum number of seepage cells', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_cellid = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'CELLID', & ! tag name
    'CELLID', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'cell identifier', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_dist = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'DIST', & ! tag name
    'DIST', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'distance', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_area = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'AREA', & ! tag name
    'AREA', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'area', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_auxvar = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'AUX', & ! tag name
    'AUXVAR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NAUX', & ! shape
    'auxiliary variables', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfspf_boundname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'BOUNDNAME', & ! tag name
    'BOUNDNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'seepage boundary name', & ! longname
    .false., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_spf_param_definitions(*) = &
    [ &
    gwfspf_auxiliary, &
    gwfspf_boundnames, &
    gwfspf_iprpak, &
    gwfspf_iprflow, &
    gwfspf_ipakcb, &
    gwfspf_some_option, &
    gwfspf_maxbound, &
    gwfspf_cellid, &
    gwfspf_dist, &
    gwfspf_area, &
    gwfspf_auxvar, &
    gwfspf_boundname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfspf_spd = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'SPF', & ! subcomponent
    'PERIOD', & ! block
    'STRESS_PERIOD_DATA', & ! tag name
    'SPD', & ! fortran variable
    'RECARRAY CELLID DIST AREA AUX BOUNDNAME', & ! type
    'MAXBOUND', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_spf_aggregate_definitions(*) = &
    [ &
    gwfspf_spd &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_spf_block_definitions(*) = &
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

end module GwfSpfInputModule

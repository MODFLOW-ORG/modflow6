! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfGncInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_gnc_param_definitions
  public gwf_gnc_aggregate_definitions
  public gwf_gnc_block_definitions
  public GwfGncParamFoundType
  public gwf_gnc_multi_package
  public gwf_gnc_subpackages

  type GwfGncParamFoundType
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: explicit = .false.
    logical :: numgnc = .false.
    logical :: numalphaj = .false.
    logical :: cellidn = .false.
    logical :: cellidm = .false.
    logical :: cellidsj = .false.
    logical :: alphasj = .false.
  end type GwfGncParamFoundType

  logical :: gwf_gnc_multi_package = .false.

  character(len=16), parameter :: &
    gwf_gnc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
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
    gwfgnc_print_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print simulated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_explicit = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'OPTIONS', & ! block
    'EXPLICIT', & ! tag name
    'EXPLICIT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use explicit GNC formulation', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_numgnc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'DIMENSIONS', & ! block
    'NUMGNC', & ! tag name
    'NUMGNC', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of ghost node corrections', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_numalphaj = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'DIMENSIONS', & ! block
    'NUMALPHAJ', & ! tag name
    'NUMALPHAJ', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of contributing factors', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_cellidn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'GNCDATA', & ! block
    'CELLIDN', & ! tag name
    'CELLIDN', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'GNC cellid n', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_cellidm = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'GNCDATA', & ! block
    'CELLIDM', & ! tag name
    'CELLIDM', & ! fortran variable
    'INTEGER1D', & ! type
    'NCELLDIM', & ! shape
    'GNC cellid n', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_cellidsj = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'GNCDATA', & ! block
    'CELLIDSJ', & ! tag name
    'CELLIDSJ', & ! fortran variable
    'INTEGER1D', & ! type
    'NUMALPHAJ', & ! shape
    'GNC contributing cells', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_alphasj = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'GNCDATA', & ! block
    'ALPHASJ', & ! tag name
    'ALPHASJ', & ! fortran variable
    'DOUBLE1D', & ! type
    'NUMALPHAJ', & ! shape
    'GNC contributing factors', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_gnc_param_definitions(*) = &
    [ &
    gwfgnc_print_input, &
    gwfgnc_print_flows, &
    gwfgnc_explicit, &
    gwfgnc_numgnc, &
    gwfgnc_numalphaj, &
    gwfgnc_cellidn, &
    gwfgnc_cellidm, &
    gwfgnc_cellidsj, &
    gwfgnc_alphasj &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfgnc_gncdata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'GNC', & ! subcomponent
    'GNCDATA', & ! block
    'GNCDATA', & ! tag name
    'GNCDATA', & ! fortran variable
    'RECARRAY CELLIDN CELLIDM CELLIDSJ ALPHASJ', & ! type
    'NUMGNC', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_gnc_aggregate_definitions(*) = &
    [ &
    gwfgnc_gncdata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_gnc_block_definitions(*) = &
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
    'GNCDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwfGncInputModule

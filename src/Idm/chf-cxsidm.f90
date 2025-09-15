! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfCxsInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_cxs_param_definitions
  public chf_cxs_aggregate_definitions
  public chf_cxs_block_definitions
  public ChfCxsParamFoundType
  public chf_cxs_multi_package
  public chf_cxs_subpackages

  type ChfCxsParamFoundType
    logical :: iprpak = .false.
    logical :: nsections = .false.
    logical :: npoints = .false.
    logical :: idcxs = .false.
    logical :: nxspoints = .false.
    logical :: xfraction = .false.
    logical :: height = .false.
    logical :: manfraction = .false.
  end type ChfCxsParamFoundType

  logical :: chf_cxs_multi_package = .false.

  character(len=16), parameter :: &
    chf_cxs_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfcxs_iprpak = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'IPRPAK', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print input to listing file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_nsections = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NSECTIONS', & ! tag name
    'NSECTIONS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of reaches', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_npoints = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NPOINTS', & ! tag name
    'NPOINTS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'total number of points defined for all reaches', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_idcxs = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'reach number for this entry', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_nxspoints = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NXSPOINTS', & ! tag name
    'NXSPOINTS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of points used to define cross section', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_xfraction = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'XFRACTION', & ! tag name
    'XFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fractional width', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_height = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'HEIGHT', & ! tag name
    'HEIGHT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'depth', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_manfraction = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'MANFRACTION', & ! tag name
    'MANFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'Mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_cxs_param_definitions(*) = &
    [ &
    chfcxs_iprpak, &
    chfcxs_nsections, &
    chfcxs_npoints, &
    chfcxs_idcxs, &
    chfcxs_nxspoints, &
    chfcxs_xfraction, &
    chfcxs_height, &
    chfcxs_manfraction &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfcxs_packagedata = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IDCXS NXSPOINTS', & ! type
    'NSECTIONS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfcxs_crosssectiondata = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'CROSSSECTIONDATA', & ! tag name
    'CROSSSECTIONDATA', & ! fortran variable
    'RECARRAY XFRACTION HEIGHT MANFRACTION', & ! type
    'NPOINTS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_cxs_aggregate_definitions(*) = &
    [ &
    chfcxs_packagedata, &
    chfcxs_crosssectiondata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_cxs_block_definitions(*) = &
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
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CROSSSECTIONDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ChfCxsInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module SwfCxsInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public swf_cxs_param_definitions
  public swf_cxs_aggregate_definitions
  public swf_cxs_block_definitions
  public SwfCxsParamFoundType
  public swf_cxs_multi_package
  public swf_cxs_subpackages

  type SwfCxsParamFoundType
    logical :: iprpak = .false.
    logical :: nsections = .false.
    logical :: npoints = .false.
    logical :: idcxs = .false.
    logical :: nxspoints = .false.
    logical :: xfraction = .false.
    logical :: height = .false.
    logical :: manfraction = .false.
  end type SwfCxsParamFoundType

  logical :: swf_cxs_multi_package = .false.

  character(len=16), parameter :: &
    swf_cxs_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfcxs_iprpak = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
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
    swfcxs_nsections = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NSECTIONS', & ! tag name
    'NSECTIONS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of reaches', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_npoints = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NPOINTS', & ! tag name
    'NPOINTS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'total number of points defined for all reaches', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_idcxs = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IDCXS', & ! tag name
    'IDCXS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'reach number for this entry', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_nxspoints = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'NXSPOINTS', & ! tag name
    'NXSPOINTS', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of points used to define cross section', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_xfraction = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'XFRACTION', & ! tag name
    'XFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fractional width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_height = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'HEIGHT', & ! tag name
    'HEIGHT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'depth', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_manfraction = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'MANFRACTION', & ! tag name
    'MANFRACTION', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'Mannings roughness coefficient', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_cxs_param_definitions(*) = &
    [ &
    swfcxs_iprpak, &
    swfcxs_nsections, &
    swfcxs_npoints, &
    swfcxs_idcxs, &
    swfcxs_nxspoints, &
    swfcxs_xfraction, &
    swfcxs_height, &
    swfcxs_manfraction &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfcxs_packagedata = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IDCXS NXSPOINTS', & ! type
    'NSECTIONS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfcxs_crosssectiondata = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'CXS', & ! subcomponent
    'CROSSSECTIONDATA', & ! block
    'CROSSSECTIONDATA', & ! tag name
    'CROSSSECTIONDATA', & ! fortran variable
    'RECARRAY XFRACTION HEIGHT MANFRACTION', & ! type
    'NPOINTS', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_cxs_aggregate_definitions(*) = &
    [ &
    swfcxs_packagedata, &
    swfcxs_crosssectiondata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    swf_cxs_block_definitions(*) = &
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

end module SwfCxsInputModule

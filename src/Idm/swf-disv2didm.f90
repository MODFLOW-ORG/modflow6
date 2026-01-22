! ** Do Not Modify! MODFLOW 6 system generated file. **
module SwfDisv2DInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public swf_disv2d_param_definitions
  public swf_disv2d_aggregate_definitions
  public swf_disv2d_block_definitions
  public SwfDisv2dParamFoundType
  public swf_disv2d_multi_package
  public swf_disv2d_subpackages

  type SwfDisv2dParamFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: grb_filerecord = .false.
    logical :: grb6 = .false.
    logical :: fileout = .false.
    logical :: grb6_filename = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: export_ascii = .false.
    logical :: crs = .false.
    logical :: nodes = .false.
    logical :: nvert = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type SwfDisv2dParamFoundType

  logical :: swf_disv2d_multi_package = .false.

  character(len=16), parameter :: &
    swf_disv2d_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_length_units = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_UNITS', & ! tag name
    'LENGTH_UNITS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model length units', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_nogrb = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'NOGRB', & ! tag name
    'NOGRB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not write binary grid file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_grb_filerecord = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB_FILERECORD', & ! tag name
    'GRB_FILERECORD', & ! fortran variable
    'RECORD GRB6 FILEOUT GRB6_FILENAME', & ! type
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
    swfdisv2d_grb6 = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6', & ! tag name
    'GRB6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'grb keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_fileout = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
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
    swfdisv2d_grb6_filename = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6_FILENAME', & ! tag name
    'GRB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of GRB information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_xorigin = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'XORIGIN', & ! tag name
    'XORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-position of the model grid origin', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_yorigin = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'YORIGIN', & ! tag name
    'YORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-position of the model grid origin', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_angrot = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'ANGROT', & ! tag name
    'ANGROT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'rotation angle', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_export_ascii = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_crs = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'OPTIONS', & ! block
    'CRS', & ! tag name
    'CRS', & ! fortran variable
    'STRING', & ! type
    'LENBIGLINE', & ! shape
    'CRS user input string', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_nodes = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NODES', & ! tag name
    'NODES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of cells per layer', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_nvert = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NVERT', & ! tag name
    'NVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_bottom = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'GRIDDATA', & ! block
    'BOTTOM', & ! tag name
    'BOTTOM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'model bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_idomain = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'GRIDDATA', & ! block
    'IDOMAIN', & ! tag name
    'IDOMAIN', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'idomain existence array', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_iv = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'VERTICES', & ! block
    'IV', & ! tag name
    'IV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'vertex number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_xv = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'VERTICES', & ! block
    'XV', & ! tag name
    'XV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_yv = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'VERTICES', & ! block
    'YV', & ! tag name
    'YV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_icell2d = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'ICELL2D', & ! tag name
    'ICELL2D', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cell2d number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_xc = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'XC', & ! tag name
    'XC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_yc = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'YC', & ! tag name
    'YC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_ncvert = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'NCVERT', & ! tag name
    'NCVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of cell vertices', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_icvert = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'ICVERT', & ! tag name
    'ICVERT', & ! fortran variable
    'INTEGER1D', & ! type
    'NCVERT', & ! shape
    'array of vertex numbers', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_disv2d_param_definitions(*) = &
    [ &
    swfdisv2d_length_units, &
    swfdisv2d_nogrb, &
    swfdisv2d_grb_filerecord, &
    swfdisv2d_grb6, &
    swfdisv2d_fileout, &
    swfdisv2d_grb6_filename, &
    swfdisv2d_xorigin, &
    swfdisv2d_yorigin, &
    swfdisv2d_angrot, &
    swfdisv2d_export_ascii, &
    swfdisv2d_crs, &
    swfdisv2d_nodes, &
    swfdisv2d_nvert, &
    swfdisv2d_bottom, &
    swfdisv2d_idomain, &
    swfdisv2d_iv, &
    swfdisv2d_xv, &
    swfdisv2d_yv, &
    swfdisv2d_icell2d, &
    swfdisv2d_xc, &
    swfdisv2d_yc, &
    swfdisv2d_ncvert, &
    swfdisv2d_icvert &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_vertices = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'VERTICES', & ! block
    'VERTICES', & ! tag name
    'VERTICES', & ! fortran variable
    'RECARRAY IV XV YV', & ! type
    'NVERT', & ! shape
    'vertices data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdisv2d_cell2d = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DISV2D', & ! subcomponent
    'CELL2D', & ! block
    'CELL2D', & ! tag name
    'CELL2D', & ! fortran variable
    'RECARRAY ICELL2D XC YC NCVERT ICVERT', & ! type
    'NODES', & ! shape
    'cell2d data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_disv2d_aggregate_definitions(*) = &
    [ &
    swfdisv2d_vertices, &
    swfdisv2d_cell2d &
    ]

  type(InputBlockDefinitionType), parameter :: &
    swf_disv2d_block_definitions(*) = &
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
    'GRIDDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'VERTICES', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CELL2D', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module SwfDisv2DInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module ChfDisv1DInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public chf_disv1d_param_definitions
  public chf_disv1d_aggregate_definitions
  public chf_disv1d_block_definitions
  public ChfDisv1dParamFoundType
  public chf_disv1d_multi_package
  public chf_disv1d_subpackages

  type ChfDisv1dParamFoundType
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
    logical :: width = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell1d = .false.
    logical :: fdc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type ChfDisv1dParamFoundType

  logical :: chf_disv1d_multi_package = .false.

  character(len=16), parameter :: &
    chf_disv1d_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_length_units = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_nogrb = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_grb_filerecord = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_grb6 = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_fileout = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_grb6_filename = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_xorigin = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'OPTIONS', & ! block
    'XORIGIN', & ! tag name
    'XORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-position origin of the model grid coordinate system', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_yorigin = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'OPTIONS', & ! block
    'YORIGIN', & ! tag name
    'YORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-position origin of the model grid coordinate system', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_angrot = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_export_ascii = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_crs = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'OPTIONS', & ! block
    'CRS', & ! tag name
    'CRS', & ! fortran variable
    'STRING', & ! type
    'LENBIGLINE', & ! shape
    'CRS user input string', & ! longname
    .false., & ! required
    .true., & ! developmode
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_nodes = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NODES', & ! tag name
    'NODES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of linear features', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_nvert = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NVERT', & ! tag name
    'NVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_width = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'GRIDDATA', & ! block
    'WIDTH', & ! tag name
    'WIDTH', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'width', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_bottom = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'GRIDDATA', & ! block
    'BOTTOM', & ! tag name
    'BOTTOM', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'bottom elevation for the one-dimensional cell', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_idomain = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_iv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_xv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_yv = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_icell1d = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'CELL1D', & ! block
    'ICELL1D', & ! tag name
    'ICELL1D', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cell1d number', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_fdc = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'CELL1D', & ! block
    'FDC', & ! tag name
    'FDC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'fractional distance to the cell center', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_ncvert = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'CELL1D', & ! block
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
    chfdisv1d_icvert = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'CELL1D', & ! block
    'ICVERT', & ! tag name
    'ICVERT', & ! fortran variable
    'INTEGER1D', & ! type
    'NCVERT', & ! shape
    'number of cell vertices', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_disv1d_param_definitions(*) = &
    [ &
    chfdisv1d_length_units, &
    chfdisv1d_nogrb, &
    chfdisv1d_grb_filerecord, &
    chfdisv1d_grb6, &
    chfdisv1d_fileout, &
    chfdisv1d_grb6_filename, &
    chfdisv1d_xorigin, &
    chfdisv1d_yorigin, &
    chfdisv1d_angrot, &
    chfdisv1d_export_ascii, &
    chfdisv1d_crs, &
    chfdisv1d_nodes, &
    chfdisv1d_nvert, &
    chfdisv1d_width, &
    chfdisv1d_bottom, &
    chfdisv1d_idomain, &
    chfdisv1d_iv, &
    chfdisv1d_xv, &
    chfdisv1d_yv, &
    chfdisv1d_icell1d, &
    chfdisv1d_fdc, &
    chfdisv1d_ncvert, &
    chfdisv1d_icvert &
    ]

  type(InputParamDefinitionType), parameter :: &
    chfdisv1d_vertices = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
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
    chfdisv1d_cell1d = InputParamDefinitionType &
    ( &
    'CHF', & ! component
    'DISV1D', & ! subcomponent
    'CELL1D', & ! block
    'CELL1D', & ! tag name
    'CELL1D', & ! fortran variable
    'RECARRAY ICELL1D FDC NCVERT ICVERT', & ! type
    'NODES', & ! shape
    'cell1d data', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    chf_disv1d_aggregate_definitions(*) = &
    [ &
    chfdisv1d_vertices, &
    chfdisv1d_cell1d &
    ]

  type(InputBlockDefinitionType), parameter :: &
    chf_disv1d_block_definitions(*) = &
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
    'CELL1D', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module ChfDisv1DInputModule

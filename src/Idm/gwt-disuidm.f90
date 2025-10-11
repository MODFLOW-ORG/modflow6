! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtDisuInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_disu_param_definitions
  public gwt_disu_aggregate_definitions
  public gwt_disu_block_definitions
  public GwtDisuParamFoundType
  public gwt_disu_multi_package
  public gwt_disu_subpackages

  type GwtDisuParamFoundType
    logical :: length_units = .false.
    logical :: nogrb = .false.
    logical :: grb_filerecord = .false.
    logical :: grb6 = .false.
    logical :: fileout = .false.
    logical :: grb6_filename = .false.
    logical :: xorigin = .false.
    logical :: yorigin = .false.
    logical :: angrot = .false.
    logical :: voffsettol = .false.
    logical :: export_ascii = .false.
    logical :: crs = .false.
    logical :: nodes = .false.
    logical :: nja = .false.
    logical :: nvert = .false.
    logical :: top = .false.
    logical :: bot = .false.
    logical :: area = .false.
    logical :: idomain = .false.
    logical :: iac = .false.
    logical :: ja = .false.
    logical :: ihc = .false.
    logical :: cl12 = .false.
    logical :: hwva = .false.
    logical :: angldegx = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type GwtDisuParamFoundType

  logical :: gwt_disu_multi_package = .false.

  character(len=16), parameter :: &
    gwt_disu_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_length_units = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_nogrb = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_grb_filerecord = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_grb6 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_fileout = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_grb6_filename = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_xorigin = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_yorigin = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_angrot = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_voffsettol = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'OPTIONS', & ! block
    'VERTICAL_OFFSET_TOLERANCE', & ! tag name
    'VOFFSETTOL', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'vertical length dimension for top and bottom checking', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_export_ascii = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_crs = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_nodes = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NODES', & ! tag name
    'NODES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of layers', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_nja = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NJA', & ! tag name
    'NJA', & ! fortran variable
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
    gwtdisu_nvert = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'DIMENSIONS', & ! block
    'NVERT', & ! tag name
    'NVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of vertices', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_top = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'TOP', & ! tag name
    'TOP', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell top elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_bot = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'BOT', & ! tag name
    'BOT', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_area = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'GRIDDATA', & ! block
    'AREA', & ! tag name
    'AREA', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'cell surface area', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_idomain = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_iac = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IAC', & ! tag name
    'IAC', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'number of cell connections', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_ja = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'JA', & ! tag name
    'JA', & ! fortran variable
    'INTEGER1D', & ! type
    'NJA', & ! shape
    'grid connectivity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_ihc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'IHC', & ! tag name
    'IHC', & ! fortran variable
    'INTEGER1D', & ! type
    'NJA', & ! shape
    'connection type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_cl12 = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'CL12', & ! tag name
    'CL12', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'connection lengths', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_hwva = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'HWVA', & ! tag name
    'HWVA', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'connection lengths', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_angldegx = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CONNECTIONDATA', & ! block
    'ANGLDEGX', & ! tag name
    'ANGLDEGX', & ! fortran variable
    'DOUBLE1D', & ! type
    'NJA', & ! shape
    'angle of face normal to connection', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_iv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_xv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_yv = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_icell2d = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_xc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_yc = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_ncvert = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwtdisu_icvert = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
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
    gwt_disu_param_definitions(*) = &
    [ &
    gwtdisu_length_units, &
    gwtdisu_nogrb, &
    gwtdisu_grb_filerecord, &
    gwtdisu_grb6, &
    gwtdisu_fileout, &
    gwtdisu_grb6_filename, &
    gwtdisu_xorigin, &
    gwtdisu_yorigin, &
    gwtdisu_angrot, &
    gwtdisu_voffsettol, &
    gwtdisu_export_ascii, &
    gwtdisu_crs, &
    gwtdisu_nodes, &
    gwtdisu_nja, &
    gwtdisu_nvert, &
    gwtdisu_top, &
    gwtdisu_bot, &
    gwtdisu_area, &
    gwtdisu_idomain, &
    gwtdisu_iac, &
    gwtdisu_ja, &
    gwtdisu_ihc, &
    gwtdisu_cl12, &
    gwtdisu_hwva, &
    gwtdisu_angldegx, &
    gwtdisu_iv, &
    gwtdisu_xv, &
    gwtdisu_yv, &
    gwtdisu_icell2d, &
    gwtdisu_xc, &
    gwtdisu_yc, &
    gwtdisu_ncvert, &
    gwtdisu_icvert &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_vertices = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'VERTICES', & ! block
    'VERTICES', & ! tag name
    'VERTICES', & ! fortran variable
    'RECARRAY IV XV YV', & ! type
    'NVERT', & ! shape
    'vertices data', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtdisu_cell2d = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'DISU', & ! subcomponent
    'CELL2D', & ! block
    'CELL2D', & ! tag name
    'CELL2D', & ! fortran variable
    'RECARRAY ICELL2D XC YC NCVERT ICVERT', & ! type
    'NODES', & ! shape
    'cell2d data', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_disu_aggregate_definitions(*) = &
    [ &
    gwtdisu_vertices, &
    gwtdisu_cell2d &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_disu_block_definitions(*) = &
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
    'CONNECTIONDATA', & ! blockname
    .true., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'VERTICES', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'CELL2D', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwtDisuInputModule

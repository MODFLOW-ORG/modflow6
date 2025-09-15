! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweDisvInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_disv_param_definitions
  public gwe_disv_aggregate_definitions
  public gwe_disv_block_definitions
  public GweDisvParamFoundType
  public gwe_disv_multi_package
  public gwe_disv_subpackages

  type GweDisvParamFoundType
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
    logical :: export_nc = .false.
    logical :: crs = .false.
    logical :: ncf_filerecord = .false.
    logical :: ncf6 = .false.
    logical :: filein = .false.
    logical :: ncf6_filename = .false.
    logical :: nlay = .false.
    logical :: ncpl = .false.
    logical :: nvert = .false.
    logical :: top = .false.
    logical :: botm = .false.
    logical :: idomain = .false.
    logical :: iv = .false.
    logical :: xv = .false.
    logical :: yv = .false.
    logical :: icell2d = .false.
    logical :: xc = .false.
    logical :: yc = .false.
    logical :: ncvert = .false.
    logical :: icvert = .false.
  end type GweDisvParamFoundType

  logical :: gwe_disv_multi_package = .false.

  character(len=16), parameter :: &
    gwe_disv_subpackages(*) = &
    [ &
    'UTL-NCF         ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwedisv_length_units = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_UNITS', & ! tag name
    'LENGTH_UNITS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model length units', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_nogrb = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'NOGRB', & ! tag name
    'NOGRB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not write binary grid file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_grb_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'GRB_FILERECORD', & ! tag name
    'GRB_FILERECORD', & ! fortran variable
    'RECORD GRB6 FILEOUT GRB6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_grb6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6', & ! tag name
    'GRB6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'grb keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_grb6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6_FILENAME', & ! tag name
    'GRB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of GRB information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_xorigin = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'XORIGIN', & ! tag name
    'XORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-position origin of the model grid coordinate system', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_yorigin = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'YORIGIN', & ! tag name
    'YORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-position origin of the model grid coordinate system', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_angrot = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'ANGROT', & ! tag name
    'ANGROT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'rotation angle', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_export_ascii = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_export_nc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_crs = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'CRS', & ! tag name
    'CRS', & ! fortran variable
    'STRING', & ! type
    'LENBIGLINE', & ! shape
    'CRS user input string', & ! longname
    .false., & ! required
    .true., & ! prerelease
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_ncf_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'NCF_FILERECORD', & ! tag name
    'NCF_FILERECORD', & ! fortran variable
    'RECORD NCF6 FILEIN NCF6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_ncf6 = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'NCF6', & ! tag name
    'NCF6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'ncf keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_ncf6_filename = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'OPTIONS', & ! block
    'NCF6_FILENAME', & ! tag name
    'NCF6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of NCF information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_nlay = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'DIMENSIONS', & ! block
    'NLAY', & ! tag name
    'NLAY', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of layers', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_ncpl = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'DIMENSIONS', & ! block
    'NCPL', & ! tag name
    'NCPL', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of cells per layer', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_nvert = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'DIMENSIONS', & ! block
    'NVERT', & ! tag name
    'NVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_top = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'GRIDDATA', & ! block
    'TOP', & ! tag name
    'TOP', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCPL', & ! shape
    'model top elevation', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_botm = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'GRIDDATA', & ! block
    'BOTM', & ! tag name
    'BOTM', & ! fortran variable
    'DOUBLE2D', & ! type
    'NCPL NLAY', & ! shape
    'model bottom elevation', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_idomain = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'GRIDDATA', & ! block
    'IDOMAIN', & ! tag name
    'IDOMAIN', & ! fortran variable
    'INTEGER2D', & ! type
    'NCPL NLAY', & ! shape
    'idomain existence array', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_iv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'VERTICES', & ! block
    'IV', & ! tag name
    'IV', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'vertex number', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_xv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'VERTICES', & ! block
    'XV', & ! tag name
    'XV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_yv = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'VERTICES', & ! block
    'YV', & ! tag name
    'YV', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for vertex', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_icell2d = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'ICELL2D', & ! tag name
    'ICELL2D', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'cell2d number', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_xc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'XC', & ! tag name
    'XC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_yc = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'YC', & ! tag name
    'YC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-coordinate for cell center', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_ncvert = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'NCVERT', & ! tag name
    'NCVERT', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of cell vertices', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_icvert = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'ICVERT', & ! tag name
    'ICVERT', & ! fortran variable
    'INTEGER1D', & ! type
    'NCVERT', & ! shape
    'array of vertex numbers', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_disv_param_definitions(*) = &
    [ &
    gwedisv_length_units, &
    gwedisv_nogrb, &
    gwedisv_grb_filerecord, &
    gwedisv_grb6, &
    gwedisv_fileout, &
    gwedisv_grb6_filename, &
    gwedisv_xorigin, &
    gwedisv_yorigin, &
    gwedisv_angrot, &
    gwedisv_export_ascii, &
    gwedisv_export_nc, &
    gwedisv_crs, &
    gwedisv_ncf_filerecord, &
    gwedisv_ncf6, &
    gwedisv_filein, &
    gwedisv_ncf6_filename, &
    gwedisv_nlay, &
    gwedisv_ncpl, &
    gwedisv_nvert, &
    gwedisv_top, &
    gwedisv_botm, &
    gwedisv_idomain, &
    gwedisv_iv, &
    gwedisv_xv, &
    gwedisv_yv, &
    gwedisv_icell2d, &
    gwedisv_xc, &
    gwedisv_yc, &
    gwedisv_ncvert, &
    gwedisv_icvert &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwedisv_vertices = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'VERTICES', & ! block
    'VERTICES', & ! tag name
    'VERTICES', & ! fortran variable
    'RECARRAY IV XV YV', & ! type
    'NVERT', & ! shape
    'vertices data', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwedisv_cell2d = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'DISV', & ! subcomponent
    'CELL2D', & ! block
    'CELL2D', & ! tag name
    'CELL2D', & ! fortran variable
    'RECARRAY ICELL2D XC YC NCVERT ICVERT', & ! type
    'NCPL', & ! shape
    'cell2d data', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_disv_aggregate_definitions(*) = &
    [ &
    gwedisv_vertices, &
    gwedisv_cell2d &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_disv_block_definitions(*) = &
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

end module GweDisvInputModule

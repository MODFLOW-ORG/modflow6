! ** Do Not Modify! MODFLOW 6 system generated file. **
module SwfDis2DInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public swf_dis2d_param_definitions
  public swf_dis2d_aggregate_definitions
  public swf_dis2d_block_definitions
  public SwfDis2dParamFoundType
  public swf_dis2d_multi_package
  public swf_dis2d_subpackages

  type SwfDis2dParamFoundType
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
    logical :: nrow = .false.
    logical :: ncol = .false.
    logical :: delr = .false.
    logical :: delc = .false.
    logical :: bottom = .false.
    logical :: idomain = .false.
  end type SwfDis2dParamFoundType

  logical :: swf_dis2d_multi_package = .false.

  character(len=16), parameter :: &
    swf_dis2d_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_length_units = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'LENGTH_UNITS', & ! tag name
    'LENGTH_UNITS', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model length units', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_nogrb = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'NOGRB', & ! tag name
    'NOGRB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'do not write binary grid file', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_grb_filerecord = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB_FILERECORD', & ! tag name
    'GRB_FILERECORD', & ! fortran variable
    'RECORD GRB6 FILEOUT GRB6_FILENAME', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_grb6 = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6', & ! tag name
    'GRB6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'grb keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_fileout = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_grb6_filename = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'GRB6_FILENAME', & ! tag name
    'GRB6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of GRB information', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_xorigin = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'XORIGIN', & ! tag name
    'XORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'x-position of the model grid origin', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_yorigin = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'YORIGIN', & ! tag name
    'YORIGIN', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'y-position of the model grid origin', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_angrot = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'ANGROT', & ! tag name
    'ANGROT', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'rotation angle', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_export_ascii = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_ASCII', & ! tag name
    'EXPORT_ASCII', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to layered ascii files.', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_crs = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'OPTIONS', & ! block
    'CRS', & ! tag name
    'CRS', & ! fortran variable
    'STRING', & ! type
    'LENBIGLINE', & ! shape
    'CRS user input string', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_nrow = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NROW', & ! tag name
    'NROW', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of rows', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_ncol = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'DIMENSIONS', & ! block
    'NCOL', & ! tag name
    'NCOL', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of columns', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_delr = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'GRIDDATA', & ! block
    'DELR', & ! tag name
    'DELR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCOL', & ! shape
    'spacing along a row', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_delc = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'GRIDDATA', & ! block
    'DELC', & ! tag name
    'DELC', & ! fortran variable
    'DOUBLE1D', & ! type
    'NROW', & ! shape
    'spacing along a column', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_bottom = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'GRIDDATA', & ! block
    'BOTTOM', & ! tag name
    'BOTTOM', & ! fortran variable
    'DOUBLE2D', & ! type
    'NCOL NROW', & ! shape
    'cell bottom elevation', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swfdis2d_idomain = InputParamDefinitionType &
    ( &
    'SWF', & ! component
    'DIS2D', & ! subcomponent
    'GRIDDATA', & ! block
    'IDOMAIN', & ! tag name
    'IDOMAIN', & ! fortran variable
    'INTEGER2D', & ! type
    'NCOL NROW', & ! shape
    'idomain existence array', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    swf_dis2d_param_definitions(*) = &
    [ &
    swfdis2d_length_units, &
    swfdis2d_nogrb, &
    swfdis2d_grb_filerecord, &
    swfdis2d_grb6, &
    swfdis2d_fileout, &
    swfdis2d_grb6_filename, &
    swfdis2d_xorigin, &
    swfdis2d_yorigin, &
    swfdis2d_angrot, &
    swfdis2d_export_ascii, &
    swfdis2d_crs, &
    swfdis2d_nrow, &
    swfdis2d_ncol, &
    swfdis2d_delr, &
    swfdis2d_delc, &
    swfdis2d_bottom, &
    swfdis2d_idomain &
    ]

  type(InputParamDefinitionType), parameter :: &
    swf_dis2d_aggregate_definitions(*) = &
    [ &
    InputParamDefinitionType &
    ( &
    '', & ! component
    '', & ! subcomponent
    '', & ! block
    '', & ! tag name
    '', & ! fortran variable
    '', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    swf_dis2d_block_definitions(*) = &
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
    ) &
    ]

end module SwfDis2DInputModule

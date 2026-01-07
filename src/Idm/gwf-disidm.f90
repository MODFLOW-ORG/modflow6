! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfDisInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_dis_param_definitions
  public gwf_dis_aggregate_definitions
  public gwf_dis_block_definitions
  public GwfDisParamFoundType
  public gwf_dis_multi_package
  public gwf_dis_subpackages

  type GwfDisParamFoundType
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
    logical :: nrow = .false.
    logical :: ncol = .false.
    logical :: delr = .false.
    logical :: delc = .false.
    logical :: top = .false.
    logical :: botm = .false.
    logical :: idomain = .false.
  end type GwfDisParamFoundType

  logical :: gwf_dis_multi_package = .false.

  character(len=16), parameter :: &
    gwf_dis_subpackages(*) = &
    [ &
    'UTL-NCF         ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfdis_length_units = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_nogrb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_grb_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_grb6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_grb6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_xorigin = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_yorigin = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_angrot = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_export_ascii = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'OPTIONS', & ! block
    'EXPORT_ARRAY_NETCDF', & ! tag name
    'EXPORT_NC', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'export array variables to netcdf output files.', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_crs = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
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
    gwfdis_ncf_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'OPTIONS', & ! block
    'NCF_FILERECORD', & ! tag name
    'NCF_FILERECORD', & ! fortran variable
    'RECORD NCF6 FILEIN NCF6_FILENAME', & ! type
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
    gwfdis_ncf6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'OPTIONS', & ! block
    'NCF6', & ! tag name
    'NCF6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'ncf keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'OPTIONS', & ! block
    'FILEIN', & ! tag name
    'FILEIN', & ! fortran variable
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
    gwfdis_ncf6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'OPTIONS', & ! block
    'NCF6_FILENAME', & ! tag name
    'NCF6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of NCF information', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_nlay = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NLAY', & ! tag name
    'NLAY', & ! fortran variable
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
    gwfdis_nrow = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NROW', & ! tag name
    'NROW', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of rows', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_ncol = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'DIMENSIONS', & ! block
    'NCOL', & ! tag name
    'NCOL', & ! fortran variable
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
    gwfdis_delr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'GRIDDATA', & ! block
    'DELR', & ! tag name
    'DELR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NCOL', & ! shape
    'spacing along a row', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_delc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'GRIDDATA', & ! block
    'DELC', & ! tag name
    'DELC', & ! fortran variable
    'DOUBLE1D', & ! type
    'NROW', & ! shape
    'spacing along a column', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_top = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'GRIDDATA', & ! block
    'TOP', & ! tag name
    'TOP', & ! fortran variable
    'DOUBLE2D', & ! type
    'NCOL NROW', & ! shape
    'cell top elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_botm = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'GRIDDATA', & ! block
    'BOTM', & ! tag name
    'BOTM', & ! fortran variable
    'DOUBLE3D', & ! type
    'NCOL NROW NLAY', & ! shape
    'cell bottom elevation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfdis_idomain = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'DIS', & ! subcomponent
    'GRIDDATA', & ! block
    'IDOMAIN', & ! tag name
    'IDOMAIN', & ! fortran variable
    'INTEGER3D', & ! type
    'NCOL NROW NLAY', & ! shape
    'idomain existence array', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_dis_param_definitions(*) = &
    [ &
    gwfdis_length_units, &
    gwfdis_nogrb, &
    gwfdis_grb_filerecord, &
    gwfdis_grb6, &
    gwfdis_fileout, &
    gwfdis_grb6_filename, &
    gwfdis_xorigin, &
    gwfdis_yorigin, &
    gwfdis_angrot, &
    gwfdis_export_ascii, &
    gwfdis_export_nc, &
    gwfdis_crs, &
    gwfdis_ncf_filerecord, &
    gwfdis_ncf6, &
    gwfdis_filein, &
    gwfdis_ncf6_filename, &
    gwfdis_nlay, &
    gwfdis_nrow, &
    gwfdis_ncol, &
    gwfdis_delr, &
    gwfdis_delc, &
    gwfdis_top, &
    gwfdis_botm, &
    gwfdis_idomain &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_dis_aggregate_definitions(*) = &
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
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_dis_block_definitions(*) = &
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

end module GwfDisInputModule

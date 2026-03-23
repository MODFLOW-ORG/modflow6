! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweNamInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_nam_param_definitions
  public gwe_nam_aggregate_definitions
  public gwe_nam_block_definitions
  public GweNamParamFoundType
  public gwe_nam_multi_package
  public gwe_nam_subpackages

  type GweNamParamFoundType
    logical :: list = .false.
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: save_flows = .false.
    logical :: idv_scale = .false.
    logical :: ncmesh2drec = .false.
    logical :: netcdf_mesh2d = .false.
    logical :: ncstructrec = .false.
    logical :: netcdf_struct = .false.
    logical :: fileout = .false.
    logical :: ncmesh2dfile = .false.
    logical :: ncstructfile = .false.
    logical :: nc_filerecord = .false.
    logical :: netcdf = .false.
    logical :: filein = .false.
    logical :: netcdf_fname = .false.
    logical :: ftype = .false.
    logical :: fname = .false.
    logical :: pname = .false.
  end type GweNamParamFoundType

  logical :: gwe_nam_multi_package = .false.

  character(len=16), parameter :: &
    gwe_nam_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwenam_list = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'LIST', & ! tag name
    'LIST', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'name of listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_print_input = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_INPUT', & ! tag name
    'PRINT_INPUT', & ! fortran variable
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
    gwenam_print_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_FLOWS', & ! tag name
    'PRINT_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'print calculated flows to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_save_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save flows for all packages to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_idv_scale = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'DEPENDENT_VARIABLE_SCALING', & ! tag name
    'IDV_SCALE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'flag to scale X and RHS', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_ncmesh2drec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NC_MESH2D_FILERECORD', & ! tag name
    'NCMESH2DREC', & ! fortran variable
    'RECORD NETCDF_MESH2D FILEOUT NCMESH2DFILE', & ! type
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
    gwenam_netcdf_mesh2d = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NETCDF_MESH2D', & ! tag name
    'NETCDF_MESH2D', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_ncstructrec = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NC_STRUCTURED_FILERECORD', & ! tag name
    'NCSTRUCTREC', & ! fortran variable
    'RECORD NETCDF_STRUCTURED FILEOUT NCSTRUCTFILE', & ! type
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
    gwenam_netcdf_struct = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NETCDF_STRUCTURED', & ! tag name
    'NETCDF_STRUCT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'budget keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_fileout = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
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
    gwenam_ncmesh2dfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NCMESH2DFILE', & ! tag name
    'NCMESH2DFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_ncstructfile = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NCSTRUCTFILE', & ! tag name
    'NCSTRUCTFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_nc_filerecord = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NC_FILERECORD', & ! tag name
    'NC_FILERECORD', & ! fortran variable
    'RECORD NETCDF FILEIN NETCDF_FILENAME', & ! type
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
    gwenam_netcdf = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NETCDF', & ! tag name
    'NETCDF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'netcdf keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
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
    gwenam_netcdf_fname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NETCDF_FILENAME', & ! tag name
    'NETCDF_FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'netcdf input filename', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_ftype = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FTYPE', & ! tag name
    'FTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'package type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_fname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'FNAME', & ! tag name
    'FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwenam_pname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PNAME', & ! tag name
    'PNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'user name for package', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_nam_param_definitions(*) = &
    [ &
    gwenam_list, &
    gwenam_print_input, &
    gwenam_print_flows, &
    gwenam_save_flows, &
    gwenam_idv_scale, &
    gwenam_ncmesh2drec, &
    gwenam_netcdf_mesh2d, &
    gwenam_ncstructrec, &
    gwenam_netcdf_struct, &
    gwenam_fileout, &
    gwenam_ncmesh2dfile, &
    gwenam_ncstructfile, &
    gwenam_nc_filerecord, &
    gwenam_netcdf, &
    gwenam_filein, &
    gwenam_netcdf_fname, &
    gwenam_ftype, &
    gwenam_fname, &
    gwenam_pname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwenam_packages = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'NAM', & ! subcomponent
    'PACKAGES', & ! block
    'PACKAGES', & ! tag name
    'PACKAGES', & ! fortran variable
    'RECARRAY FTYPE FNAME PNAME', & ! type
    '', & ! shape
    'package list', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_nam_aggregate_definitions(*) = &
    [ &
    gwenam_packages &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_nam_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGES', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GweNamInputModule

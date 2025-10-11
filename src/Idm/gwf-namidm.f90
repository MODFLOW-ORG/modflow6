! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfNamInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_nam_param_definitions
  public gwf_nam_aggregate_definitions
  public gwf_nam_block_definitions
  public GwfNamParamFoundType
  public gwf_nam_multi_package
  public gwf_nam_subpackages

  type GwfNamParamFoundType
    logical :: list = .false.
    logical :: print_input = .false.
    logical :: print_flows = .false.
    logical :: save_flows = .false.
    logical :: newtonoptions = .false.
    logical :: newton = .false.
    logical :: under_relaxation = .false.
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
  end type GwfNamParamFoundType

  logical :: gwf_nam_multi_package = .false.

  character(len=16), parameter :: &
    gwf_nam_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfnam_list = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_print_input = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_print_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_save_flows = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_newtonoptions = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTONOPTIONS', & ! tag name
    'NEWTONOPTIONS', & ! fortran variable
    'RECORD NEWTON UNDER_RELAXATION', & ! type
    '', & ! shape
    'newton keyword and options', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnam_newton = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'NEWTON', & ! tag name
    'NEWTON', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson formulation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnam_under_relaxation = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'NAM', & ! subcomponent
    'OPTIONS', & ! block
    'UNDER_RELAXATION', & ! tag name
    'UNDER_RELAXATION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to activate Newton-Raphson UNDER_RELAXATION option', & ! longname
    .false., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfnam_ncmesh2drec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_netcdf_mesh2d = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_ncstructrec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_netcdf_struct = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_ncmesh2dfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_ncstructfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_nc_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_netcdf = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_netcdf_fname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_ftype = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_fname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwfnam_pname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwf_nam_param_definitions(*) = &
    [ &
    gwfnam_list, &
    gwfnam_print_input, &
    gwfnam_print_flows, &
    gwfnam_save_flows, &
    gwfnam_newtonoptions, &
    gwfnam_newton, &
    gwfnam_under_relaxation, &
    gwfnam_ncmesh2drec, &
    gwfnam_netcdf_mesh2d, &
    gwfnam_ncstructrec, &
    gwfnam_netcdf_struct, &
    gwfnam_fileout, &
    gwfnam_ncmesh2dfile, &
    gwfnam_ncstructfile, &
    gwfnam_nc_filerecord, &
    gwfnam_netcdf, &
    gwfnam_filein, &
    gwfnam_netcdf_fname, &
    gwfnam_ftype, &
    gwfnam_fname, &
    gwfnam_pname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfnam_packages = InputParamDefinitionType &
    ( &
    'GWF', & ! component
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
    gwf_nam_aggregate_definitions(*) = &
    [ &
    gwfnam_packages &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_nam_block_definitions(*) = &
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

end module GwfNamInputModule

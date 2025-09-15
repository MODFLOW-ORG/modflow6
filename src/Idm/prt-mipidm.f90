! ** Do Not Modify! MODFLOW 6 system generated file. **
module PrtMipInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public prt_mip_param_definitions
  public prt_mip_aggregate_definitions
  public prt_mip_block_definitions
  public PrtMipParamFoundType
  public prt_mip_multi_package
  public prt_mip_subpackages

  type PrtMipParamFoundType
    logical :: export_ascii = .false.
    logical :: porosity = .false.
    logical :: retfactor = .false.
    logical :: izone = .false.
  end type PrtMipParamFoundType

  logical :: prt_mip_multi_package = .false.

  character(len=16), parameter :: &
    prt_mip_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtmip_export_ascii = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'MIP', & ! subcomponent
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
    prtmip_porosity = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'MIP', & ! subcomponent
    'GRIDDATA', & ! block
    'POROSITY', & ! tag name
    'POROSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'porosity', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtmip_retfactor = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'MIP', & ! subcomponent
    'GRIDDATA', & ! block
    'RETFACTOR', & ! tag name
    'RETFACTOR', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'retardation factor', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtmip_izone = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'MIP', & ! subcomponent
    'GRIDDATA', & ! block
    'IZONE', & ! tag name
    'IZONE', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'zone number', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prt_mip_param_definitions(*) = &
    [ &
    prtmip_export_ascii, &
    prtmip_porosity, &
    prtmip_retfactor, &
    prtmip_izone &
    ]

  type(InputParamDefinitionType), parameter :: &
    prt_mip_aggregate_definitions(*) = &
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
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    ) &
    ]

  type(InputBlockDefinitionType), parameter :: &
    prt_mip_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
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

end module PrtMipInputModule

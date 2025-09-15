! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfStoInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_sto_param_definitions
  public gwf_sto_aggregate_definitions
  public gwf_sto_block_definitions
  public GwfStoParamFoundType
  public gwf_sto_multi_package
  public gwf_sto_subpackages

  type GwfStoParamFoundType
    logical :: ipakcb = .false.
    logical :: istor_coef = .false.
    logical :: ss_confined_only = .false.
    logical :: tvs_filerecord = .false.
    logical :: tvs6 = .false.
    logical :: filein = .false.
    logical :: tvs6_filename = .false.
    logical :: export_ascii = .false.
    logical :: export_nc = .false.
    logical :: iorig_ss = .false.
    logical :: iconf_ss = .false.
    logical :: iconvert = .false.
    logical :: ss = .false.
    logical :: sy = .false.
    logical :: steady_state = .false.
    logical :: transient = .false.
  end type GwfStoParamFoundType

  logical :: gwf_sto_multi_package = .false.

  character(len=16), parameter :: &
    gwf_sto_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfsto_ipakcb = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'IPAKCB', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to save NPF flows', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_istor_coef = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'STORAGECOEFFICIENT', & ! tag name
    'ISTOR_COEF', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate SS is read as storage coefficient', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_ss_confined_only = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'SS_CONFINED_ONLY', & ! tag name
    'SS_CONFINED_ONLY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'keyword to indicate specific storage only applied under confined&
     & conditions', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_tvs_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'TVS_FILERECORD', & ! tag name
    'TVS_FILERECORD', & ! fortran variable
    'RECORD TVS6 FILEIN TVS6_FILENAME', & ! type
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
    gwfsto_tvs6 = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'TVS6', & ! tag name
    'TVS6', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'tvs keyword', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_filein = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
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
    gwfsto_tvs6_filename = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'TVS6_FILENAME', & ! tag name
    'TVS6_FILENAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name of TVS information', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_export_ascii = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
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
    gwfsto_export_nc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
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
    gwfsto_iorig_ss = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_ORIGINAL_SPECIFIC_STORAGE', & ! tag name
    'IORIG_SS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'development option for original specific storage', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_iconf_ss = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_OLDSTORAGEFORMULATION', & ! tag name
    'ICONF_SS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'development option flag for old storage formulation', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_iconvert = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'GRIDDATA', & ! block
    'ICONVERT', & ! tag name
    'ICONVERT', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'convertible indicator', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_ss = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'GRIDDATA', & ! block
    'SS', & ! tag name
    'SS', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'specific storage', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_sy = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'GRIDDATA', & ! block
    'SY', & ! tag name
    'SY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'specific yield', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_steady_state = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'PERIOD', & ! block
    'STEADY-STATE', & ! tag name
    'STEADY_STATE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'steady state indicator', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfsto_transient = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'STO', & ! subcomponent
    'PERIOD', & ! block
    'TRANSIENT', & ! tag name
    'TRANSIENT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'transient indicator', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_sto_param_definitions(*) = &
    [ &
    gwfsto_ipakcb, &
    gwfsto_istor_coef, &
    gwfsto_ss_confined_only, &
    gwfsto_tvs_filerecord, &
    gwfsto_tvs6, &
    gwfsto_filein, &
    gwfsto_tvs6_filename, &
    gwfsto_export_ascii, &
    gwfsto_export_nc, &
    gwfsto_iorig_ss, &
    gwfsto_iconf_ss, &
    gwfsto_iconvert, &
    gwfsto_ss, &
    gwfsto_sy, &
    gwfsto_steady_state, &
    gwfsto_transient &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_sto_aggregate_definitions(*) = &
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
    gwf_sto_block_definitions(*) = &
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
    ), &
    InputBlockDefinitionType( &
    'PERIOD', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .true. & ! block_variable
    ) &
    ]

end module GwfStoInputModule

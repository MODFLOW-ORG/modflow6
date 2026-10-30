! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfVsfInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_vsf_param_definitions
  public gwf_vsf_aggregate_definitions
  public gwf_vsf_block_definitions
  public GwfVsfParamFoundType
  public gwf_vsf_multi_package
  public gwf_vsf_is_advanced
  public gwf_vsf_subpackages

  type GwfVsfParamFoundType
    logical :: soil_model = .false.
    logical :: dev_soil_param = .false.
    logical :: model_kr = .false.
    logical :: storage_scheme = .false.
    logical :: kr_averaging = .false.
    logical :: sat_filerecord = .false.
    logical :: saturation = .false.
    logical :: fileout = .false.
    logical :: saturationfile = .false.
    logical :: phead_filerecord = .false.
    logical :: pressurehead = .false.
    logical :: pheadfile = .false.
    logical :: iunsat = .false.
    logical :: porosity = .false.
    logical :: satres = .false.
    logical :: alphahvk = .false.
    logical :: nhvk = .false.
    logical :: betahvk = .false.
    logical :: khvk = .false.
    logical :: alphavgn = .false.
    logical :: nvgn = .false.
  end type GwfVsfParamFoundType

  logical :: gwf_vsf_multi_package = .false.
  logical :: gwf_vsf_is_advanced = .false.

  character(len=16), parameter :: &
    gwf_vsf_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_soil_model = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'SOIL_MODEL', & ! tag name
    'SOIL_MODEL', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model used for the soil characteristic functions', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_dev_soil_param = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_SOIL_PARAM', & ! tag name
    'DEV_SOIL_PARAM', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'set a soil parameterization', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_model_kr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'MODEL_KR', & ! tag name
    'MODEL_KR', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'model used for the relative permeability', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_storage_scheme = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'STORAGE_SCHEME', & ! tag name
    'STORAGE_SCHEME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'discretization method for storage calculation', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_kr_averaging = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'KR_AVERAGING', & ! tag name
    'KR_AVERAGING', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'relative permeability averaging method', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_sat_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'SAT_FILERECORD', & ! tag name
    'SAT_FILERECORD', & ! fortran variable
    'RECORD SATURATION FILEOUT SATURATIONFILE', & ! type
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
    gwfvsf_saturation = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'SATURATION', & ! tag name
    'SATURATION', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'saturation keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'FILEOUT', & ! tag name
    'FILEOUT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'fileout keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_saturationfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'SATURATIONFILE', & ! tag name
    'SATURATIONFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'saturation file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_phead_filerecord = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'PHEAD_FILERECORD', & ! tag name
    'PHEAD_FILERECORD', & ! fortran variable
    'RECORD PRESSUREHEAD FILEOUT PHEADFILE', & ! type
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
    gwfvsf_pressurehead = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'PRESSUREHEAD', & ! tag name
    'PRESSUREHEAD', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'pressure head keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_pheadfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'OPTIONS', & ! block
    'PHEADFILE', & ! tag name
    'PHEADFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'pressure head file keyword', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_iunsat = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'IUNSAT', & ! tag name
    'IUNSAT', & ! fortran variable
    'INTEGER1D', & ! type
    'NODES', & ! shape
    'unsaturated zone indicator', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_porosity = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'POROSITY', & ! tag name
    'POROSITY', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'porosity', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_satres = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'SATRES', & ! tag name
    'SATRES', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'residual saturation', & ! longname
    .true., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_alphahvk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'ALPHAHVK', & ! tag name
    'ALPHAHVK', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'alpha', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_nhvk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'NHVK', & ! tag name
    'NHVK', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'n', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_betahvk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'BETAHVK', & ! tag name
    'BETAHVK', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'beta', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_khvk = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'KHVK', & ! tag name
    'KHVK', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'k', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_alphavgn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'ALPHAVGN', & ! tag name
    'ALPHAVGN', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'alpha', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfvsf_nvgn = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'VSF', & ! subcomponent
    'GRIDDATA', & ! block
    'NVGN', & ! tag name
    'NVGN', & ! fortran variable
    'DOUBLE1D', & ! type
    'NODES', & ! shape
    'n', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .true., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_vsf_param_definitions(*) = &
    [ &
    gwfvsf_soil_model, &
    gwfvsf_dev_soil_param, &
    gwfvsf_model_kr, &
    gwfvsf_storage_scheme, &
    gwfvsf_kr_averaging, &
    gwfvsf_sat_filerecord, &
    gwfvsf_saturation, &
    gwfvsf_fileout, &
    gwfvsf_saturationfile, &
    gwfvsf_phead_filerecord, &
    gwfvsf_pressurehead, &
    gwfvsf_pheadfile, &
    gwfvsf_iunsat, &
    gwfvsf_porosity, &
    gwfvsf_satres, &
    gwfvsf_alphahvk, &
    gwfvsf_nhvk, &
    gwfvsf_betahvk, &
    gwfvsf_khvk, &
    gwfvsf_alphavgn, &
    gwfvsf_nvgn &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwf_vsf_aggregate_definitions(*) = &
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
    gwf_vsf_block_definitions(*) = &
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

end module GwfVsfInputModule

! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwfBuyInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwf_buy_param_definitions
  public gwf_buy_aggregate_definitions
  public gwf_buy_block_definitions
  public GwfBuyParamFoundType
  public gwf_buy_multi_package
  public gwf_buy_subpackages

  type GwfBuyParamFoundType
    logical :: hhform_rhs = .false.
    logical :: denseref = .false.
    logical :: density_fr = .false.
    logical :: density = .false.
    logical :: fileout = .false.
    logical :: densityfile = .false.
    logical :: dev_efh_form = .false.
    logical :: nrhospecies = .false.
    logical :: irhospec = .false.
    logical :: drhodc = .false.
    logical :: crhoref = .false.
    logical :: modelname = .false.
    logical :: auxspeciesname = .false.
  end type GwfBuyParamFoundType

  logical :: gwf_buy_multi_package = .false.

  character(len=16), parameter :: &
    gwf_buy_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_hhform_rhs = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'HHFORMULATION_RHS', & ! tag name
    'HHFORM_RHS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'hh formulation on right-hand side', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_denseref = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'DENSEREF', & ! tag name
    'DENSEREF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reference density', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_density_fr = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'DENSITY_FILERECORD', & ! tag name
    'DENSITY_FR', & ! fortran variable
    'RECORD DENSITY FILEOUT DENSITYFILE', & ! type
    '', & ! shape
    '', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_density = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'DENSITY', & ! tag name
    'DENSITY', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'density keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_fileout = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
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
    gwfbuy_densityfile = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'DENSITYFILE', & ! tag name
    'DENSITYFILE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file keyword', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_dev_efh_form = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_EFH_FORMULATION', & ! tag name
    'DEV_EFH_FORM', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'use equivalent freshwater head formulation', & ! longname
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_nrhospecies = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'DIMENSIONS', & ! block
    'NRHOSPECIES', & ! tag name
    'NRHOSPECIES', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'number of species used in density equation of state', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_irhospec = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'IRHOSPEC', & ! tag name
    'IRHOSPEC', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    'species number for this entry', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_drhodc = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'DRHODC', & ! tag name
    'DRHODC', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'slope of the density-concentration line', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_crhoref = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'CRHOREF', & ! tag name
    'CRHOREF', & ! fortran variable
    'DOUBLE', & ! type
    '', & ! shape
    'reference concentration value', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_modelname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'MODELNAME', & ! tag name
    'MODELNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'modelname', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_auxspeciesname = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'AUXSPECIESNAME', & ! tag name
    'AUXSPECIESNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'auxspeciesname', & ! longname
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_buy_param_definitions(*) = &
    [ &
    gwfbuy_hhform_rhs, &
    gwfbuy_denseref, &
    gwfbuy_density_fr, &
    gwfbuy_density, &
    gwfbuy_fileout, &
    gwfbuy_densityfile, &
    gwfbuy_dev_efh_form, &
    gwfbuy_nrhospecies, &
    gwfbuy_irhospec, &
    gwfbuy_drhodc, &
    gwfbuy_crhoref, &
    gwfbuy_modelname, &
    gwfbuy_auxspeciesname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwfbuy_packagedata = InputParamDefinitionType &
    ( &
    'GWF', & ! component
    'BUY', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY IRHOSPEC DRHODC CRHOREF MODELNAME AUXSPECIESNAME', & ! type
    'NRHOSPECIES', & ! shape
    '', & ! longname
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwf_buy_aggregate_definitions(*) = &
    [ &
    gwfbuy_packagedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwf_buy_block_definitions(*) = &
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
    'PACKAGEDATA', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GwfBuyInputModule

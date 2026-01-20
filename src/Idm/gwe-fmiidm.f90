! ** Do Not Modify! MODFLOW 6 system generated file. **
module GweFmiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwe_fmi_param_definitions
  public gwe_fmi_aggregate_definitions
  public gwe_fmi_block_definitions
  public GweFmiParamFoundType
  public gwe_fmi_multi_package
  public gwe_fmi_subpackages

  type GweFmiParamFoundType
    logical :: save_flows = .false.
    logical :: imbalancecorrect = .false.
    logical :: flowtype = .false.
    logical :: filein = .false.
    logical :: fname = .false.
  end type GweFmiParamFoundType

  logical :: gwe_fmi_multi_package = .false.

  character(len=16), parameter :: &
    gwe_fmi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwefmi_save_flows = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flow imbalance correction to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwefmi_imbalancecorrect = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_IMBALANCE_CORRECTION', & ! tag name
    'IMBALANCECORRECT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'correct for flow imbalance', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwefmi_flowtype = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
    'FLOWTYPE', & ! tag name
    'FLOWTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'flow type', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwefmi_filein = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
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
    gwefmi_fname = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
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
    gwe_fmi_param_definitions(*) = &
    [ &
    gwefmi_save_flows, &
    gwefmi_imbalancecorrect, &
    gwefmi_flowtype, &
    gwefmi_filein, &
    gwefmi_fname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwefmi_packagedata = InputParamDefinitionType &
    ( &
    'GWE', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY FLOWTYPE FILEIN FNAME', & ! type
    '', & ! shape
    'flowtype list', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwe_fmi_aggregate_definitions(*) = &
    [ &
    gwefmi_packagedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwe_fmi_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PACKAGEDATA', & ! blockname
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module GweFmiInputModule

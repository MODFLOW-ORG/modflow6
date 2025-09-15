! ** Do Not Modify! MODFLOW 6 system generated file. **
module GwtFmiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public gwt_fmi_param_definitions
  public gwt_fmi_aggregate_definitions
  public gwt_fmi_block_definitions
  public GwtFmiParamFoundType
  public gwt_fmi_multi_package
  public gwt_fmi_subpackages

  type GwtFmiParamFoundType
    logical :: save_flows = .false.
    logical :: imbalancecorrect = .false.
    logical :: flowtype = .false.
    logical :: filein = .false.
    logical :: fname = .false.
  end type GwtFmiParamFoundType

  logical :: gwt_fmi_multi_package = .false.

  character(len=16), parameter :: &
    gwt_fmi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtfmi_save_flows = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save calculated flow imbalance correction to budget file', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtfmi_imbalancecorrect = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'OPTIONS', & ! block
    'FLOW_IMBALANCE_CORRECTION', & ! tag name
    'IMBALANCECORRECT', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'correct for flow imbalance', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtfmi_flowtype = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
    'FLOWTYPE', & ! tag name
    'FLOWTYPE', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'flow type', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwtfmi_filein = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
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
    gwtfmi_fname = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
    'FNAME', & ! tag name
    'FNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    'file name', & ! longname
    .true., & ! required
    .false., & ! prerelease
    .true., & ! multi-record
    .true., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_fmi_param_definitions(*) = &
    [ &
    gwtfmi_save_flows, &
    gwtfmi_imbalancecorrect, &
    gwtfmi_flowtype, &
    gwtfmi_filein, &
    gwtfmi_fname &
    ]

  type(InputParamDefinitionType), parameter :: &
    gwtfmi_packagedata = InputParamDefinitionType &
    ( &
    'GWT', & ! component
    'FMI', & ! subcomponent
    'PACKAGEDATA', & ! block
    'PACKAGEDATA', & ! tag name
    'PACKAGEDATA', & ! fortran variable
    'RECARRAY FLOWTYPE FILEIN FNAME', & ! type
    '', & ! shape
    'flowtype list', & ! longname
    .false., & ! required
    .false., & ! prerelease
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    gwt_fmi_aggregate_definitions(*) = &
    [ &
    gwtfmi_packagedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    gwt_fmi_block_definitions(*) = &
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

end module GwtFmiInputModule

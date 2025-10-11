! ** Do Not Modify! MODFLOW 6 system generated file. **
module PrtFmiInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public prt_fmi_param_definitions
  public prt_fmi_aggregate_definitions
  public prt_fmi_block_definitions
  public PrtFmiParamFoundType
  public prt_fmi_multi_package
  public prt_fmi_subpackages

  type PrtFmiParamFoundType
    logical :: save_flows = .false.
    logical :: flowtype = .false.
    logical :: filein = .false.
    logical :: fname = .false.
  end type PrtFmiParamFoundType

  logical :: prt_fmi_multi_package = .false.

  character(len=16), parameter :: &
    prt_fmi_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtfmi_save_flows = InputParamDefinitionType &
    ( &
    'PRT', & ! component
    'FMI', & ! subcomponent
    'OPTIONS', & ! block
    'SAVE_FLOWS', & ! tag name
    'SAVE_FLOWS', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'save cell-by-cell flows to budget file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    prtfmi_flowtype = InputParamDefinitionType &
    ( &
    'PRT', & ! component
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
    prtfmi_filein = InputParamDefinitionType &
    ( &
    'PRT', & ! component
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
    prtfmi_fname = InputParamDefinitionType &
    ( &
    'PRT', & ! component
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
    prt_fmi_param_definitions(*) = &
    [ &
    prtfmi_save_flows, &
    prtfmi_flowtype, &
    prtfmi_filein, &
    prtfmi_fname &
    ]

  type(InputParamDefinitionType), parameter :: &
    prtfmi_packagedata = InputParamDefinitionType &
    ( &
    'PRT', & ! component
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
    prt_fmi_aggregate_definitions(*) = &
    [ &
    prtfmi_packagedata &
    ]

  type(InputBlockDefinitionType), parameter :: &
    prt_fmi_block_definitions(*) = &
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

end module PrtFmiInputModule

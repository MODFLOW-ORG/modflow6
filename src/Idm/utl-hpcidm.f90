! ** Do Not Modify! MODFLOW 6 system generated file. **
module UtlHpcInputModule
  use ConstantsModule, only: LENVARNAME
  use InputDefinitionModule, only: InputParamDefinitionType, &
                                   InputBlockDefinitionType
  private
  public utl_hpc_param_definitions
  public utl_hpc_aggregate_definitions
  public utl_hpc_block_definitions
  public UtlHpcParamFoundType
  public utl_hpc_multi_package

  type UtlHpcParamFoundType
    logical :: dev_log_mpi = .false.
    logical :: mname = .false.
    logical :: mrank = .false.
  end type UtlHpcParamFoundType

  logical :: utl_hpc_multi_package = .false.

  type(InputParamDefinitionType), parameter :: &
    utlhpc_dev_log_mpi = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'HPC', & ! subcomponent
    'OPTIONS', & ! block
    'DEV_LOG_MPI', & ! tag name
    'DEV_LOG_MPI', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    .false., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlhpc_mname = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'HPC', & ! subcomponent
    'PARTITIONS', & ! block
    'MNAME', & ! tag name
    'MNAME', & ! fortran variable
    'STRING', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utlhpc_mrank = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'HPC', & ! subcomponent
    'PARTITIONS', & ! block
    'MRANK', & ! tag name
    'MRANK', & ! fortran variable
    'INTEGER', & ! type
    '', & ! shape
    .true., & ! required
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_hpc_param_definitions(*) = &
    [ &
    utlhpc_dev_log_mpi, &
    utlhpc_mname, &
    utlhpc_mrank &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlhpc_partitions = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'HPC', & ! subcomponent
    'PARTITIONS', & ! block
    'PARTITIONS', & ! tag name
    'PARTITIONS', & ! fortran variable
    'RECARRAY MNAME MRANK', & ! type
    '', & ! shape
    .true., & ! required
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_hpc_aggregate_definitions(*) = &
    [ &
    utlhpc_partitions &
    ]

  type(InputBlockDefinitionType), parameter :: &
    utl_hpc_block_definitions(*) = &
    [ &
    InputBlockDefinitionType( &
    'OPTIONS', & ! blockname
    .false., & ! required
    .false., & ! aggregate
    .false. & ! block_variable
    ), &
    InputBlockDefinitionType( &
    'PARTITIONS', & ! blockname
    .true., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module UtlHpcInputModule

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
  public utl_hpc_subpackages

  type UtlHpcParamFoundType
    logical :: print_table = .false.
    logical :: dev_log_mpi = .false.
    logical :: mname = .false.
    logical :: mrank = .false.
  end type UtlHpcParamFoundType

  logical :: utl_hpc_multi_package = .false.

  character(len=16), parameter :: &
    utl_hpc_subpackages(*) = &
    [ &
    '                ' &
    ]

  type(InputParamDefinitionType), parameter :: &
    utlhpc_print_table = InputParamDefinitionType &
    ( &
    'UTL', & ! component
    'HPC', & ! subcomponent
    'OPTIONS', & ! block
    'PRINT_TABLE', & ! tag name
    'PRINT_TABLE', & ! fortran variable
    'KEYWORD', & ! type
    '', & ! shape
    'model print table to listing file', & ! longname
    .false., & ! required
    .false., & ! developmode
    .false., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

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
    'log mpi traffic', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    'model name', & ! longname
    .true., & ! required
    .false., & ! developmode
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
    'model rank', & ! longname
    .true., & ! required
    .false., & ! developmode
    .true., & ! multi-record
    .false., & ! preserve case
    .false., & ! layered
    .false. & ! timeseries
    )

  type(InputParamDefinitionType), parameter :: &
    utl_hpc_param_definitions(*) = &
    [ &
    utlhpc_print_table, &
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
    'list of partition numbers', & ! longname
    .false., & ! required
    .false., & ! developmode
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
    .false., & ! required
    .true., & ! aggregate
    .false. & ! block_variable
    ) &
    ]

end module UtlHpcInputModule

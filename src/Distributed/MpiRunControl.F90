module MpiRunControlModule
#if defined(__WITH_PETSC__)
#include <petsc/finclude/petscksp.h>
  ! cannot use all symbols because of clash with 'use mpi'
  use petscksp, only: PETSC_COMM_WORLD, PetscInitialize, PetscFinalize, &
                      PETSC_NULL_CHARACTER, PETSC_NULL_OPTIONS
#endif
  use mpi
  use MpiWorldModule
  use SimVariablesModule, only: proc_id, nr_procs
  use SimStagesModule
  use KindModule, only: I4B, LGP, DP
  use STLVecIntModule
  use NumericalSolutionModule
  use RunControlModule, only: RunControlType
  implicit none
  private

  public :: create_mpi_run_control

  type, public, extends(RunControlType) :: MpiRunControlType
  contains
    ! override
    procedure :: start => mpi_ctrl_start
    procedure :: finish => mpi_ctrl_finish
    procedure :: after_con_cr => mpi_ctrl_after_con_cr
    ! private
    procedure, private :: wait_for_debugger
  end type MpiRunControlType

contains

  function create_mpi_run_control() result(controller)
    class(RunControlType), pointer :: controller
    ! local
    class(MpiRunControlType), pointer :: mpi_controller

    allocate (mpi_controller)
    controller => mpi_controller

  end function create_mpi_run_control

  subroutine mpi_ctrl_start(this)
    use ErrorUtilModule, only: pstop_alternative
    class(MpiRunControlType) :: this
    ! local
    integer :: ierr
    character(len=*), parameter :: petsc_db_file = '.petscrc'
    logical(LGP) :: petsc_db_exists, wait_dbg, is_parallel_mode
    type(MpiWorldType), pointer :: mpi_world

    ! set mpi abort function
    pstop_alternative => mpi_stop

    wait_dbg = .false.
    mpi_world => get_mpi_world()

    ! if PETSc we need their initialize
#if defined(__WITH_PETSC__)
    ! PetscInitialize calls MPI_Init only when it is not called yet,
    ! which could be through the API. If it is already called, we
    ! should assign the MPI communicator to PETSC_COMM_WORLD first
    ! (PETSc manual)
    if (mpi_world%has_comm()) then
      PETSC_COMM_WORLD = mpi_world%comm
    end if

    inquire (file=petsc_db_file, exist=petsc_db_exists)
    if (.not. petsc_db_exists) then
      call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
      CHKERRQ(ierr)
    else
      call PetscInitialize(petsc_db_file, ierr)
      CHKERRQ(ierr)
    end if

    if (.not. mpi_world%has_comm()) then
      call mpi_world%set_comm(PETSC_COMM_WORLD)
    end if

    call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-wait_dbg', wait_dbg, ierr)
    CHKERRQ(ierr)
    call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
                             '-p', is_parallel_mode, ierr)
    CHKERRQ(ierr)
#else
    if (.not. mpi_world%has_comm()) then
      call MPI_Init(ierr)
      call CHECK_MPI(ierr)
      call mpi_world%set_comm(MPI_COMM_WORLD)
    end if
#endif

    call mpi_world%init()

    call MPI_Comm_size(mpi_world%comm, nr_procs, ierr)
    call MPI_Comm_rank(mpi_world%comm, proc_id, ierr)

    ! possibly wait to attach debugger here
    if (wait_dbg) call this%wait_for_debugger()

    ! start everything else by calling parent
    call this%RunControlType%start()

  end subroutine mpi_ctrl_start

  subroutine wait_for_debugger(this)
    class(MpiRunControlType) :: this
    ! local
    integer :: ierr
    integer(I4B) :: icnt
    type(MpiWorldType), pointer :: mpi_world

    mpi_world => get_mpi_world()
    if (proc_id == 0) then
      icnt = 0
      write (*, *) 'Hit enter to continue...'
      read (*, *)
    end if
    call MPI_Barrier(mpi_world%comm, ierr)

  end subroutine wait_for_debugger

  subroutine mpi_ctrl_finish(this)
    use ErrorUtilModule, only: pstop_alternative
    class(MpiRunControlType) :: this
    ! local
    integer :: ierr

    ! release MPI related memory in router before MPI_Finalize
    call this%virtual_data_mgr%router%finalize()

    ! finish mpi
#if defined(__WITH_PETSC__)
    ! NB: PetscFinalize calls MPI_Finalize only when MPI_Init
    ! was called before PetscInitialize
    call PetscFinalize(ierr)
    CHKERRQ(ierr)
#else
    call MPI_Finalize(ierr)
    call CHECK_MPI(ierr)
#endif

    pstop_alternative => null()

    ! finish everything else by calling parent
    call this%RunControlType%finish()

  end subroutine mpi_ctrl_finish

  !> @brief Actions after creating connections
  !<
  subroutine mpi_ctrl_after_con_cr(this)
    use VirtualDataListsModule
    use VirtualModelModule, only: VirtualModelType, &
                                  get_virtual_model_from_list, &
                                  get_virtual_model
    use VirtualExchangeModule, only: VirtualExchangeType, &
                                     get_virtual_exchange_from_list, &
                                     get_virtual_exchange
    class(MpiRunControlType) :: this
    ! local
    integer(I4B) :: i, j, id, irank
    integer(I4B) :: nr_models, nr_exgs, nr_remotes, max_nr_remotes
    type(STLVecInt) :: remote_models, remote_exgs
    integer(I4B), dimension(:, :), pointer :: remote_models_per_process
    integer(I4B), dimension(:, :), pointer :: remote_exgs_per_process
    class(VirtualModelType), pointer :: vm
    class(VirtualExchangeType), pointer :: ve
    type(MpiWorldType), pointer :: mpi_world
    integer :: ierr

    mpi_world => get_mpi_world()

    ! activate halo through base
    call this%RunControlType%after_con_cr()

    ! compose list of remote models/exchanges to receive
    call remote_models%init()
    nr_models = virtual_model_list%Count()
    do i = 1, nr_models
      vm => get_virtual_model_from_list(virtual_model_list, i)
      if (vm%is_active .and. .not. vm%is_local) then
        ! remote and active
        call remote_models%push_back(vm%id)
      end if
    end do
    call remote_exgs%init()
    nr_exgs = virtual_exchange_list%Count()
    do i = 1, nr_exgs
      ve => get_virtual_exchange_from_list(virtual_exchange_list, i)
      if (ve%is_active .and. .not. ve%is_local) then
        ! remote and active
        call remote_exgs%push_back(ve%id)
      end if
    end do

    ! Models: find max for allocation
    nr_remotes = remote_models%size
    call MPI_Allreduce(nr_remotes, max_nr_remotes, 1, MPI_INTEGER, MPI_MAX, &
                       mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    allocate (remote_models_per_process(max_nr_remotes, nr_procs))
    remote_models_per_process = 0

    ! Models: fill local portion and reduce
    do i = 1, remote_models%size
      remote_models_per_process(i, proc_id + 1) = remote_models%at(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, remote_models_per_process, &
                       max_nr_remotes * nr_procs, MPI_INTEGER, MPI_MAX, &
                       mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    ! Models: set remotes to virtual models
    do i = 1, nr_procs
      do j = 1, max_nr_remotes
        id = remote_models_per_process(j, i)
        if (id > 0) then
          ! assign zero-based rank number to virtual model
          vm => get_virtual_model(id)
          if (vm%is_local) then
            ! only for local models
            irank = i - 1
            call vm%rcv_ranks%push_back_unique(irank)
          end if
        end if
      end do
    end do

    ! Exchanges: find max for allocation
    nr_remotes = remote_exgs%size
    call MPI_Allreduce(nr_remotes, max_nr_remotes, 1, MPI_INTEGER, MPI_MAX, &
                       mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    allocate (remote_exgs_per_process(max_nr_remotes, nr_procs))
    remote_exgs_per_process = 0

    ! Exchanges: fill local portion and reduce
    do i = 1, remote_exgs%size
      remote_exgs_per_process(i, proc_id + 1) = remote_exgs%at(i)
    end do
    call MPI_Allreduce(MPI_IN_PLACE, remote_exgs_per_process, &
                       max_nr_remotes * nr_procs, MPI_INTEGER, MPI_MAX, &
                       mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    ! Exchanges: set remotes to virtual exchanges
    do i = 1, nr_procs
      do j = 1, max_nr_remotes
        id = remote_exgs_per_process(j, i)
        if (id > 0) then
          ! assign zero-based rank number to virtual exchange
          ve => get_virtual_exchange(id)
          if (ve%is_local) then
            ! only for local exchanges
            irank = i - 1
            call ve%rcv_ranks%push_back_unique(irank)
          end if
        end if
      end do
    end do

    ! clean up
    call remote_models%destroy()
    call remote_exgs%destroy()

    deallocate (remote_models_per_process)
    deallocate (remote_exgs_per_process)

  end subroutine mpi_ctrl_after_con_cr

end module MpiRunControlModule

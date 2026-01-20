module RunControlModule
  use KindModule, only: I4B
  use SimStagesModule
  use VirtualDataManagerModule
  use MapperModule
  use ListsModule, only: basesolutionlist
  use NumericalSolutionModule, only: NumericalSolutionType
  use ProfilerModule
  implicit none
  private

  public :: create_seq_run_control

  type, public :: RunControlType
    class(VirtualDataManagerType), pointer :: virtual_data_mgr !< syncs globally accessible data, timely, by
                                                               !! linking (local) or message passing (remote)
    type(MapperType) :: mapper !< a 'mapper' for copying data between two memory addresses
  contains
    procedure :: start => ctrl_start
    procedure :: at_stage => ctrl_at_stage
    procedure :: finish => ctrl_finish
    procedure :: after_con_cr => ctrl_after_con_cr

    ! private
    procedure, private :: init_handler
    procedure, private :: before_con_df
    procedure, private :: after_con_df
    procedure, private :: destroy
  end type RunControlType

contains

  function create_seq_run_control() result(run_controller)
    class(RunControlType), pointer :: run_controller

    allocate (run_controller)

  end function create_seq_run_control

  subroutine ctrl_start(this)
    class(RunControlType) :: this

    allocate (this%virtual_data_mgr)

  end subroutine ctrl_start

  subroutine ctrl_finish(this)
    use SimVariablesModule, only: iout
    use MemoryManagerModule, only: mem_write_usage, mem_da
    use TimerModule, only: elapsed_time
    use SimModule, only: final_message
    class(RunControlType) :: this

    ! clean up
    call this%destroy()

    ! -- Write memory usage, elapsed time and terminate
    call mem_write_usage(iout)
    call mem_da()

    ! stop and print timings
    call g_prof%stop(g_prof%tmr_finalize)
    call g_prof%stop(g_prof%tmr_run)
    call g_prof%print(iout)
    call g_prof%destroy()

    call elapsed_time(iout, 1)
    call final_message()

  end subroutine ctrl_finish

  !> @brief This will call the handler for a particular stage
  !< in the simulation run
  subroutine ctrl_at_stage(this, stage)
    class(RunControlType) :: this
    integer(I4B) :: stage

    if (stage == STG_BFR_MDL_DF) then
      call this%init_handler()
    else if (stage == STG_AFT_CON_CR) then
      call this%after_con_cr()
    else if (stage == STG_BFR_CON_DF) then
      call this%before_con_df()
    else if (stage == STG_AFT_CON_DF) then
      call this%after_con_df()
    end if

    call this%virtual_data_mgr%synchronize(stage)
    call this%mapper%scatter(0, stage)

  end subroutine ctrl_at_stage

  subroutine init_handler(this)
    use SimVariablesModule, only: simulation_mode
    class(RunControlType), target :: this

    call this%virtual_data_mgr%create(simulation_mode)
    call this%virtual_data_mgr%init()
    call this%mapper%init()

  end subroutine init_handler

  !> @brief Actions after connections have been created
  !<
  subroutine ctrl_after_con_cr(this)
    class(RunControlType) :: this

    call this%virtual_data_mgr%activate_halo()

  end subroutine ctrl_after_con_cr

  !> @brief Actions before defining the connections
  !!
  !! Set up the virtual data manager:
  !! The models and exchanges in the halo for this interface
  !! have been determined. Add them to the virtual data manager
  !! for synchronization. (After which the interface model
  !< grids can be constructed)
  subroutine before_con_df(this)
    class(RunControlType), target :: this
    ! local
    integer(I4B) :: i
    class(*), pointer :: obj_ptr
    class(NumericalSolutionType), pointer :: sol

    ! Add (halo) models and exchanges to the virtual
    ! solutions. Set the synchronization handler
    ! in the numerical solution.
    do i = 1, basesolutionlist%Count()
      obj_ptr => basesolutionlist%GetItem(i)
      select type (obj_ptr)
      class is (NumericalSolutionType)
        sol => obj_ptr
        call this%virtual_data_mgr%add_solution(sol)
        sol%synchronize => rc_solution_sync
        sol%synchronize_ctx => this
      end select
    end do

    ! The remote data fields in exchanges need to
    ! be copied in from the virtual exchanges
    call this%mapper%add_exchange_vars()

  end subroutine before_con_df

  !> @brief Actions after defining connections
  !<
  subroutine after_con_df(this)
    class(RunControlType) :: this

    ! Reduce the halo
    call this%virtual_data_mgr%compress_halo()

    ! Add variables in interface models to the mapper
    call this%mapper%add_interface_vars()

  end subroutine after_con_df

  !> @brief Synchronizes from within numerical solution (delegate)
  !<
  subroutine rc_solution_sync(num_sol, stage, ctx)
    use NumericalSolutionModule, only: NumericalSolutionType
    class(NumericalSolutionType) :: num_sol
    integer(I4B) :: stage
    class(*), pointer :: ctx

    select type (ctx)
    class is (RunControlType)
      call ctx%virtual_data_mgr%synchronize_sln(num_sol%id, stage)
      call ctx%mapper%scatter(num_sol%id, stage)
    end select

  end subroutine rc_solution_sync

  subroutine destroy(this)
    class(RunControlType) :: this

    call this%virtual_data_mgr%destroy()
    deallocate (this%virtual_data_mgr)

  end subroutine destroy

end module RunControlModule

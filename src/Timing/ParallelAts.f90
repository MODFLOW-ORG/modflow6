module ParallelAtsModule
  use mpi
  use SimVariablesModule, only: proc_id
  use MpiWorldModule, only: get_mpi_world, MpiWorldType, CHECK_MPI
  use KindModule, only: I4B, DP
  use AdaptiveTimeStepModule
  implicit none
  private

  public :: ParallelAtsType

  !> @brief Extends ATS so we can do synchronization
  !< of time step size (including dtstable!) across ranks
  type, extends(AtsType) :: ParallelAtsType
  contains
    procedure, public :: ats_submit_delt => par_ats_submit_delt
    procedure, public :: ats_set_delt => par_ats_set_delt
  end type ParallelAtsType

contains

  !> @ brief Allow and external caller to submit preferred time step
  !!
  !!  Submit a preferred time step length.  Alternatively, if idir is
  !!  is passed, then either increase or decrease the submitted time
  !!  step by the dtadj input variable.
  !!
  !<
  subroutine par_ats_submit_delt(this, kstp, kper, dt, sloc, idir)
    class(ParallelAtsType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    real(DP), intent(in) :: dt
    character(len=*), intent(in) :: sloc
    integer(I4B), intent(in), optional :: idir
    ! local
    type(MpiWorldType), pointer :: mpi_world
    integer :: ierr
    real(DP) :: global_dt

    if (.not. this%isAdaptivePeriod(kper)) return

    ! reduce dtstable over all ranks
    mpi_world => get_mpi_world()
    call MPI_Allreduce(this%dtstable, global_dt, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    this%dtstable = global_dt

    call this%AtsType%ats_submit_delt(kstp, kper, dt, sloc, idir)

  end subroutine par_ats_submit_delt

  !> @ brief Set time step
  !!
  !! Set the time step length (delt) for this time step using the ATS
  !! controls, but with synchronization across ranks.
  !<
  subroutine par_ats_set_delt(this, kstp, kper, pertim, perlencurrent, delt)
    class(ParallelAtsType) :: this
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: kper
    real(DP), intent(inout) :: pertim
    real(DP), intent(in) :: perlencurrent
    real(DP), intent(inout) :: delt
    ! local
    type(MpiWorldType), pointer :: mpi_world
    integer :: ierr
    real(DP) :: global_delt

    ! reduce delt over all ranks
    mpi_world => get_mpi_world()
    call MPI_Allreduce(delt, global_delt, 1, MPI_DOUBLE_PRECISION, &
                       MPI_MIN, mpi_world%comm, ierr)
    call CHECK_MPI(ierr)

    delt = global_delt

    call this%AtsType%ats_set_delt(kstp, kper, pertim, perlencurrent, delt)

  end subroutine par_ats_set_delt

end module ParallelAtsModule

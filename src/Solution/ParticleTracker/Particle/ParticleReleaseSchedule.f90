!> @brief Particle release scheduling.
module ParticleReleaseScheduleModule

  use ArrayHandlersModule, only: ExpandArray
  use ConstantsModule, only: DZERO, DONE, LINELENGTH
  use KindModule, only: I4B, LGP, DP
  use MathUtilModule, only: is_close
  use TimeSelectModule, only: TimeSelectType
  use TimeStepSelectModule, only: TimeStepSelectType

  implicit none
  private
  public :: ParticleReleaseScheduleType
  public :: create_release_schedule

  !> @brief Particle release scheduling utility.
  !!
  !! The release schedule composes a time selection object for any
  !! explicitly specified release times, and a time step selection
  !! object for release times specified in period/time step terms.
  !!
  !! Release time coincidence is computed within a given tolerance;
  !! times closer than the tolerance are merged into a single time.
  !!
  !! The release schedule must be refreshed each time step. This is
  !! achieved by calling `advance()`. After this, the `times` member
  !! is a debounced/consolidated schedule for the current time step.
  !<
  type :: ParticleReleaseScheduleType
    real(DP), allocatable :: times(:) !< release times
    real(DP) :: tolerance !< release time coincidence tolerance
    type(TimeSelectType), pointer :: time_select !< time selection
    type(TimeStepSelectType), pointer :: step_select !< time step selection
  contains
    procedure :: advance
    procedure :: any
    procedure :: count
    procedure :: deallocate
    procedure :: log
    procedure :: schedule
  end type ParticleReleaseScheduleType

contains

  !> @brief Create a new release schedule object.
  function create_release_schedule(tol) result(sched)
    real(DP), intent(in) :: tol !< coincident release time tolerance
    type(ParticleReleaseScheduleType), pointer :: sched !< schedule pointer

    allocate (sched)
    allocate (sched%times(0))
    allocate (sched%time_select)
    allocate (sched%step_select)
    call sched%time_select%init()
    call sched%step_select%init()
    sched%tolerance = tol

  end function create_release_schedule

  !> @brief Deallocate the release schedule.
  subroutine deallocate (this)
    class(ParticleReleaseScheduleType), intent(inout) :: this !< this instance

    deallocate (this%times)
    call this%time_select%deallocate()
    call this%step_select%deallocate()
    deallocate (this%time_select)
    deallocate (this%step_select)

  end subroutine deallocate

  !> @brief Write the release schedule to the given output unit.
  subroutine log(this, iout)
    class(ParticleReleaseScheduleType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: iout !< output unit
    character(len=*), parameter :: fmt = &
      &"(6x,A,': ',50(G0,' '))"

    if (this%any()) then
      write (iout, fmt) 'RELEASE SCHEDULE', this%times
    else
      write (iout, "(1x,a,1x,a)") 'NO RELEASES SCHEDULED'
    end if
  end subroutine log

  !> @brief Add a release time to the schedule.
  !!
  !! To schedule multiple release times at once, expand
  !! and populate the time selection object by hand. DO
  !! NOT attempt to manipulate the times array; this is
  !! a read-only property which the schedule maintains.
  !<
  subroutine schedule(this, trelease)
    class(ParticleReleaseScheduleType), intent(inout) :: this
    real(DP), intent(in) :: trelease
    call ExpandArray(this%times)
    this%times(size(this%times)) = trelease
  end subroutine schedule

  !> @brief Refresh the schedule for the current time step.
  !!
  !! This involves several tasks: first, advance the time
  !! selection. Then, if period-block release setting lines
  !! are provided, reinitialize the time step selection for
  !! the given period. Finally, refresh the schedule array,
  !! deduplicating any times closer than the set tolerance.
  !!
  !! This routine is idempotent.
  !<
  subroutine advance(this, lines)
    use TdisModule, only: totimc, kstp, endofperiod
    class(ParticleReleaseScheduleType), intent(inout) :: this
    character(len=LINELENGTH), intent(in), optional :: lines(:)
    integer(I4B) :: it, i
    real(DP) :: tprevious
    real(DP) :: trelease

    ! Advance the time selection.
    call this%time_select%advance()

    ! Reinitialize the time step selection if new
    ! period-block release settings are provided.
    if (present(lines)) then
      call this%step_select%init()
      do i = 1, size(lines)
        call this%step_select%read(lines(i))
      end do
    end if

    ! Reinitialize the release time schedule.
    if (allocated(this%times)) deallocate (this%times)
    allocate (this%times(0))

    tprevious = -DONE
    trelease = -DONE

    ! Add a release time configured by period-block
    ! settings, if one is scheduled this time step.
    if (this%step_select%is_selected(kstp, endofperiod=endofperiod)) then
      trelease = totimc
      call this%schedule(trelease)
      tprevious = trelease
    end if

    ! Schedule explicitly specified release times, up
    ! to the configured tolerance of coincidence
    if (this%time_select%any()) then
      do it = this%time_select%selection(1), this%time_select%selection(2)
        trelease = this%time_select%times(it)
        if (tprevious >= DZERO .and. is_close( &
            tprevious, &
            trelease, &
            atol=this%tolerance)) cycle

        call this%schedule(trelease)
        tprevious = trelease
      end do
    end if
  end subroutine advance

  !> @brief Check if any releases are scheduled.
  !!
  !! Note: be sure to call advance() before calling this function,
  !! or the result may still be associated with a prior time step.
  !<
  logical function any(this) result(a)
    class(ParticleReleaseScheduleType) :: this
    a = this%count() > 0
  end function any

  !> @brief Return the number of releases scheduled.
  !!
  !! Note: be sure to call advance() before calling this function,
  !! or the result may still be associated with a prior time step.
  !<
  integer function count(this) result(n)
    class(ParticleReleaseScheduleType) :: this
    n = size(this%times)
  end function count

end module ParticleReleaseScheduleModule

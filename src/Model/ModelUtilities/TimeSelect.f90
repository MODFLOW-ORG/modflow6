!> @brief Specify times for some event to occur.
module TimeSelectModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE
  use ArrayHandlersModule, only: ExpandArray
  use ErrorUtilModule, only: pstop
  use SortModule, only: qsort

  implicit none
  public :: TimeSelectType

  !> @brief Represents a series of instants at which some event should occur.
  !!
  !! Maintains an array of configured times which can be sliced to match e.g.
  !! the current period & time step. Slicing can be performed manually, with
  !! the select() routine, or automatically, with the advance() routine, for
  !! a convenient view onto the applicable subset of the complete time array.
  !!
  !! Time selection uses the interval convention (t0, t1] (exclusive lower
  !! bound, inclusive upper bound). This ensures times at exact time step
  !! boundaries are captured by the earlier time step without duplication.
  !!
  !! Array storage can be expanded manually. Note: array expansion must take
  !! place before selection; when expand() is called the selection is wiped.
  !! Alternatively, the extend() routine will automatically expand the array
  !! and sort it.
  !!
  !! Most use cases likely assume a strictly increasing time selection; this
  !! can be checked with increasing(). Note that the sort() routine does not
  !! check for duplicates, and should usually be followed by an increasing()
  !! check before the time selection is used.
  !<
  type :: TimeSelectType
    real(DP), allocatable :: times(:)
    integer(I4B) :: selection(2)
  contains
    procedure :: deallocate
    procedure :: expand
    procedure :: init
    procedure :: increasing
    procedure :: log
    procedure :: select
    procedure :: advance
    procedure :: any
    procedure :: count
    procedure :: sort
    procedure :: extend
  end type TimeSelectType

contains

  !> @brief Deallocate the time selection object.
  subroutine deallocate (this)
    class(TimeSelectType) :: this
    deallocate (this%times)
  end subroutine deallocate

  !> @brief Expand capacity by the given amount. Resets the current slice.
  subroutine expand(this, increment)
    class(TimeSelectType) :: this
    integer(I4B), optional, intent(in) :: increment

    call ExpandArray(this%times, increment=increment)
    this%selection = (/1, size(this%times)/)
  end subroutine expand

  !> @brief Initialize or clear the time selection object.
  subroutine init(this)
    class(TimeSelectType) :: this

    if (allocated(this%times)) deallocate (this%times)
    allocate (this%times(0))
    this%selection = (/0, 0/)
  end subroutine

  !> @brief Determine if times strictly increase.
  !!
  !! Returns true if the times array strictly increases,
  !! as well as if the times array is empty, or not yet
  !! allocated. Note that this function operates on the
  !! entire times array, not the current selection. Note
  !! also that this function conducts exact comparisons;
  !! deduplication with tolerance must be done manually.
  !<
  function increasing(this) result(inc)
    class(TimeSelectType) :: this
    logical(LGP) :: inc
    integer(I4B) :: i
    real(DP) :: l, t

    inc = .true.
    if (.not. allocated(this%times)) return
    do i = 1, size(this%times)
      t = this%times(i)
      if (i /= 1) then
        if (l >= t) then
          inc = .false.
          return
        end if
      end if
      l = t
    end do
  end function increasing

  !> @brief Show the current time selection, if any.
  subroutine log(this, iout, verb)
    ! dummy
    class(TimeSelectType) :: this !< this instance
    integer(I4B), intent(in) :: iout !< output unit
    character(len=*), intent(in) :: verb !< selection name
    ! formats
    character(len=*), parameter :: fmt = &
      &"(6x,'THE FOLLOWING TIMES WILL BE ',A,': ',50(G0,' '))"

    if (this%any()) then
      write (iout, fmt) verb, this%times(this%selection(1):this%selection(2))
    else
      write (iout, "(a,1x,a)") 'NO TIMES WILL BE', verb
    end if
  end subroutine log

  !> @brief Select times in the interval (t0, t1] (exclusive lower, inclusive upper).
  !!
  !! Finds and stores the index of the first time strictly after the start time,
  !! and of the last time at or before the end time. This interval convention
  !! ensures that times at exact time step boundaries are captured by the later
  !! (earlier in simulation time) step without duplication.
  !!
  !! Allows filtering the times for e.g. a particular stress period and time step.
  !! Array indices are assumed to start at 1. If no times are found to fall within
  !! the selection (i.e. the interval falls entirely between two consecutive times
  !! or beyond the time range), indices are set to [-1, -1].
  !!
  !! The given start and end times are first checked against currently stored
  !! indices to avoid recalculating them if possible, allowing multiple consuming
  !! components (e.g., subdomain particle tracking solutions) to share the object
  !! efficiently, provided all proceed through stress periods and time steps in
  !! lockstep, i.e. they all solve any given period/step before any will proceed
  !! to the next.
  !<
  subroutine select(this, t0, t1, changed)
    ! dummy
    class(TimeSelectType) :: this
    real(DP), intent(in) :: t0, t1
    logical(LGP), intent(inout), optional :: changed
    ! local
    integer(I4B) :: i, i0, i1
    integer(I4B) :: l, u, lp, up
    real(DP) :: t

    ! by default, need to iterate over all times
    i0 = 1
    i1 = size(this%times)

    ! if no times fall within the slice, set to [-1, -1]
    l = -1
    u = -1

    ! previous bounding indices
    lp = this%selection(1)
    up = this%selection(2)

    ! Check if we can reuse either the lower or upper bound.
    ! The lower doesn't need to change if it indexes the first
    ! time strictly after the slice's start (i.e., the previous
    ! time is at or before t0, and this time is after t0).
    ! The upper doesn't need to change if it indexes the last
    ! time at or before the slice's end.
    if (lp > 0 .and. up > 0) then
      if (lp > 1) then
        if (this%times(lp - 1) <= t0 .and. &
            this%times(lp) > t0) then
          l = lp
          i0 = l
        end if
      end if
      if (up > 1 .and. up < i1) then
        if (this%times(up + 1) > t1 .and. &
            this%times(up) <= t1) then
          u = up
          i1 = u
        end if
      end if
      if (l == lp .and. u == up) then
        this%selection = (/l, u/)
        if (present(changed)) changed = .false.
        return
      end if
    end if

    ! recompute bounding indices if needed
    do i = i0, i1
      t = this%times(i)
      if (l < 0 .and. t > t0 .and. t <= t1) l = i
      if (l > 0 .and. t <= t1) u = i
    end do
    this%selection = (/l, u/)
    if (present(changed)) changed = l /= lp .or. u /= up

  end subroutine

  !> @brief Update the selection to the current time step.
  subroutine advance(this)
    ! modules
    use TdisModule, only: kper, kstp, nper, nstp, totimc, delt
    ! dummy
    class(TimeSelectType) :: this
    ! local
    real(DP) :: l, u

    if (kper == 1 .and. kstp == 1) then
      ! For first time step, use a small negative lower bound
      ! capture times at t=0.0 despite exclusive lower bound
      l = -epsilon(DZERO)
    else
      l = totimc
    end if
    if (kper == nper .and. kstp == nstp(kper)) then
      ! For last time step, use a large upper bound to
      ! capture times beyond the end of the simulation
      u = huge(DONE)
    else
      u = totimc + delt
    end if
    call this%select(l, u)
  end subroutine advance

  !> @brief Check if any times are currently selected.
  !!
  !! Indicates whether any times are selected for the
  !! current time step.
  !!
  !! Note that this routine does NOT indicate whether
  !! the times array has nonzero size; use the size
  !! intrinsic for that.
  !<
  function any(this) result(a)
    class(TimeSelectType) :: this
    logical(LGP) :: a

    a = all(this%selection > 0)
  end function any

  !> @brief Return the number of times currently selected.
  !!
  !! Returns the number of times selected for the current
  !! time step.
  !!
  !! Note that this routine does NOT return the total size
  !! of the times array; use the size intrinsic for that.
  !<
  function count(this) result(n)
    class(TimeSelectType) :: this
    integer(I4B) :: n

    if (this%any()) then
      n = this%selection(2) - this%selection(1)
    else
      n = 0
    end if
  end function count

  !> @brief Sort the time selection in increasing order.
  !!
  !! Note that this routine does NOT remove duplicate times.
  !! Call increasing() to check for duplicates in the array.
  !<
  subroutine sort(this)
    class(TimeSelectType) :: this
    integer(I4B), allocatable :: indx(:)

    allocate (indx(size(this%times)))
    call qsort(indx, this%times)
    deallocate (indx)
  end subroutine sort

  !> @brief Extend the time selection with the given array.
  !!
  !! This routine sorts the selection after appending the
  !! elements of the given array, but users should likely
  !! still call increasing() to check for duplicate times.
  !<
  subroutine extend(this, a)
    class(TimeSelectType) :: this
    real(DP) :: a(:)

    this%times = [this%times, a]
    call this%sort()
  end subroutine extend

end module TimeSelectModule

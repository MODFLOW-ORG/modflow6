!> @brief Time step selection module.
module TimeStepSelectModule

  use KindModule, only: DP, I4B, LGP
  use ArrayHandlersModule, only: expandarray
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use InputOutputModule, only: urword

  implicit none
  private
  public :: TimeStepSelectType

  !> @brief Time step selection type.
  !!
  !! Represents a selection of time steps as configured in an input file's
  !! period block settings. The object should be initiated with the init()
  !! procedure. The read() procedure accepts a character string of form:
  !!
  !!   ALL
  !!   STEPS 1 4 5 6
  !!   FIRST
  !!   LAST
  !!   FREQUENCY 4
  !!
  !! The read() procedure may be invoked multiple times to select multiple
  !! time steps. Note that a character string reusing a keyword which has
  !! been used for a previous read() invocation will override the previous
  !! setting using that keyword. To combine multiple settings, be sure the
  !! keywords are different on each invocation, e.g.:
  !!
  !!   FIRST
  !!   LAST
  !!   STEPS 2
  !!
  !! The is_selected() function indicates whether the given time step is
  !! active. This function accepts an optional argument, indicating that
  !! the time step is the last in the stress period.
  !<
  type :: TimeStepSelectType
    logical(LGP) :: all
    logical(LGP) :: first
    logical(LGP) :: last
    integer(I4B) :: freq
    integer(I4B), allocatable :: steps(:)
  contains
    procedure :: deallocate
    procedure :: init
    procedure :: log
    procedure :: read
    procedure :: is_selected
    procedure :: any
  end type TimeStepSelectType

contains

  !> @brief Deallocate the time step selection object.
  subroutine deallocate (this)
    class(TimeSTepSelectType) :: this
    deallocate (this%steps)
  end subroutine deallocate

  !> @brief Initialize the time step selection object.
  subroutine init(this)
    class(TimeStepSelectType) :: this !< this instance

    if (allocated(this%steps)) deallocate (this%steps)
    allocate (this%steps(0))
    this%freq = 0
    this%first = .false.
    this%last = .false.
    this%all = .false.
  end subroutine init

  subroutine log(this, iout, verb)
    ! dummy
    class(TimeStepSelectType) :: this !< this instance
    integer(I4B), intent(in) :: iout !< output unit
    character(len=*), intent(in) :: verb !< selection name
    ! formats
    character(len=*), parameter :: fmt_steps = &
      &"(6x,'THE FOLLOWING STEPS WILL BE ',A,': ',50(I0,' '))"
    character(len=*), parameter :: fmt_freq = &
      &"(6x,'THE FOLLOWING FREQUENCY WILL BE ',A,': ',I0)"

    if (this%all) then
      write (iout, "(6x,a,a)") 'ALL TIME STEPS WILL BE ', verb
    end if
    if (size(this%steps) > 0) then
      write (iout, fmt_steps) verb, this%steps
    end if
    if (this%freq > 0) then
      write (iout, fmt_freq) verb, this%freq
    end if
    if (this%first) then
      write (iout, "(6x,a,a)") 'THE FIRST TIME STEP WILL BE ', verb
    end if
    if (this%last) then
      write (iout, "(6x,a,a)") 'THE LAST TIME STEP WILL BE ', verb
    end if
  end subroutine log

  !> @brief Read a line of input and prepare the selection object.
  subroutine read (this, line)
    class(TimeStepSelectType) :: this !< this instance
    character(len=*), intent(in) :: line !< input line

    character(len=len(line)) :: l
    integer(I4B) :: n, lloc, istart, istop, ival
    real(DP) :: rval

    l(:) = line(:)
    lloc = 1

    call urword(l, lloc, istart, istop, 1, ival, rval, 0, 0)
    select case (l(istart:istop))
    case ('ALL')
      this%all = .true.
    case ('STEPS')
      listsearch: do
        call urword(l, lloc, istart, istop, 2, ival, rval, -1, 0)
        if (ival > 0) then
          n = size(this%steps)
          call expandarray(this%steps)
          this%steps(n + 1) = ival
          cycle listsearch
        end if
        exit listsearch
      end do listsearch
    case ('FREQUENCY')
      call urword(l, lloc, istart, istop, 2, ival, rval, -1, 0)
      this%freq = ival
    case ('FIRST')
      this%first = .true.
    case ('LAST')
      this%last = .true.
    case default
      write (errmsg, '(2a)') &
        'Looking for ALL, STEPS, FIRST, LAST, OR FREQUENCY. Found: ', &
        trim(adjustl(line))
      call store_error(errmsg, terminate=.TRUE.)
    end select
  end subroutine read

  !> @brief Indicates whether the given time step is selected.
  logical function is_selected(this, kstp, endofperiod)
    ! dummy
    class(TimeStepSelectType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in), optional :: endofperiod !< whether last step of stress period
    ! local
    integer(I4B) :: i, n
    logical(LGP) :: lend

    if (present(endofperiod)) then
      lend = endofperiod
    else
      lend = .false.
    end if

    is_selected = .false.
    if (this%all) is_selected = .true.
    if (kstp == 1 .and. this%first) is_selected = .true.
    if (lend .and. this%last) is_selected = .true.
    if (this%freq > 0) then
      if (mod(kstp, this%freq) == 0) is_selected = .true.
    end if
    n = size(this%steps)
    if (n > 0) then
      do i = 1, n
        if (kstp == this%steps(i)) then
          is_selected = .true.
          exit
        end if
      end do
    end if
  end function is_selected

  !> @brief Indicates whether any time steps are selected.
  logical function any(this) result(a)
    class(TimeStepSelectType) :: this !< this instance
    a = (this%all .or. &
         this%first .or. &
         this%last .or. &
         this%freq > 0 .or. &
         size(this%steps) > 0)
  end function any

end module TimeStepSelectModule

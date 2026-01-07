!> @brief Print/save manager module.
module PrintSaveManagerModule

  use KindModule, only: DP, I4B, LGP
  use ArrayHandlersModule, only: expandarray
  use SimModule, only: store_error
  use InputOutputModule, only: urword
  use TimeStepSelectModule, only: TimeStepSelectType

  implicit none
  private
  public :: PrintSaveManagerType
  public :: create_psm

  !> @brief Print/save manager type.
  !!
  !! Stores user settings as configured in an input file's period block
  !! and determines whether data should be printed to a list (log) file
  !! or saved to disk.
  !!
  !! The object should be initiated with the init() procedure.
  !!
  !! The rp() procedure will read a character string and configure the
  !! manager, where the character string may be of the following form:
  !!
  !!   PRINT ALL
  !!   PRINT STEPS 1 4 5 6
  !!   PRINT FIRST
  !!   PRINT LAST
  !!   PRINT FREQUENCY 4
  !!   SAVE ALL
  !!   SAVE STEPS 1 4 5 6
  !!   SAVE FIRST
  !!   SAVE LAST
  !!   SAVE FREQUENCY 4
  !!
  !! The should_print() and should_save() functions indicate whether
  !! to save or print during the current time step.
  !<
  type :: PrintSaveManagerType
    private
    type(TimeStepSelectType), pointer, public :: save_steps
    type(TimeStepSelectType), pointer, public :: print_steps
  contains
    procedure :: allocate
    procedure :: deallocate
    procedure :: init
    procedure :: read
    procedure :: set
    procedure :: should_print
    procedure :: should_save
  end type PrintSaveManagerType

contains

  !> @brief Initialize or clear the print/save manager.
  function create_psm() result(psm)
    type(PrintSaveManagerType), pointer :: psm !< the print/save manager
    allocate (psm)
    call psm%allocate()
  end function create_psm

  subroutine allocate (this)
    class(PrintSaveManagerType) :: this
    allocate (this%save_steps)
    allocate (this%print_steps)
    call this%save_steps%init()
    call this%print_steps%init()
  end subroutine allocate

  subroutine deallocate (this)
    class(PrintSaveManagerType) :: this !< this instance
    if (associated(this%save_steps)) deallocate (this%save_steps)
    if (associated(this%print_steps)) deallocate (this%print_steps)
  end subroutine deallocate

  subroutine init(this)
    class(PrintSaveManagerType) :: this
    call this%deallocate()
    call this%allocate()
  end subroutine init

  !> @ brief Read a line of input and prepare the manager.
  subroutine read (this, linein, iout)
    ! dummy
    class(PrintSaveManagerType) :: this !< this instance
    character(len=*), intent(in) :: linein !< input line
    integer(I4B), intent(in) :: iout !< output file unit
    ! local
    character(len=len(linein)) :: line
    integer(I4B) :: lloc, istart, istop, ival
    real(DP) :: rval

    line(:) = linein(:)
    lloc = 1
    call urword(line, lloc, istart, istop, 1, ival, rval, 0, 0)

    select case (line(istart:istop))
    case ('PRINT')
      call this%print_steps%read(line(istop + 2:))
      if (iout > 0) &
        call this%print_steps%log(iout, verb="PRINTED")
    case ('SAVE')
      call this%save_steps%read(line(istop + 2:))
      if (iout > 0) &
        call this%save_steps%log(iout, verb="SAVED")
    case default
      call store_error( &
        'Looking for PRINT or SAVE. Found: '//trim(adjustl(line)), &
        terminate=.TRUE.)
    end select
  end subroutine read

  !> @ brief set data from an input record.
  subroutine set(this, ocaction, ocsetting, iout)
    ! dummy
    class(PrintSaveManagerType) :: this !< this instance
    character(len=*) :: ocaction
    character(len=*) :: ocsetting
    integer(I4B), intent(in) :: iout !< output file unit
    select case (ocaction)
    case ('PRINT')
      call this%print_steps%read(ocsetting)
      if (iout > 0) &
        call this%print_steps%log(iout, verb="PRINTED")
    case ('SAVE')
      call this%save_steps%read(ocsetting)
      if (iout > 0) &
        call this%save_steps%log(iout, verb="SAVED")
    case default
      call store_error( &
        'Looking for PRINT or SAVE. Found: '//trim(adjustl(ocaction)), &
        terminate=.TRUE.)
    end select
  end subroutine set

  !> @ brief Determine if printing is enabled for this time step.
  logical function should_print(this, kstp, endofperiod)
    class(PrintSaveManagerType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< whether last step of stress period

    should_print = this%print_steps%is_selected(kstp, endofperiod)
  end function should_print

  !> @ brief Determine if saving is enabled for this time step.
  logical function should_save(this, kstp, endofperiod)
    class(PrintSaveManagerType) :: this !< this instance
    integer(I4B), intent(in) :: kstp !< current time step
    logical(LGP), intent(in) :: endofperiod !< whether last step of stress period

    should_save = this%save_steps%is_selected(kstp, endofperiod)
  end function should_save

end module PrintSaveManagerModule

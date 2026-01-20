!> @brief This module contains the LongLineReaderType
!!
!! The LongLineReader is a utility for reading text lines
!! from mf6 input files.  It calls u9rdcom (which calls
!! get_line) to read the first non-commented line of an
!! input file.  The LongLineReader can emulate the Fortran
!! backspace command by calling the bkspc method, which stores
!! the current line in last_line, and will return last_line
!! upon the next call to rdcom.  The LongLineReader was
!! implemented to replace all Fortran backspace calls, due
!! to a bug in ifort and ifx that prevented the backspace
!! command from working properly with non-advancing IO.
!!
!<
module LongLineReaderModule

  use, intrinsic :: iso_fortran_env, only: IOSTAT_END
  use KindModule, only: I4B
  use SimModule, only: store_error
  use InputOutputModule, only: u9rdcom

  implicit none

  private
  public :: LongLineReaderType

  !> @brief LongLineReaderType
  !!
  !!  Object for reading input from mf6 input files
  !!
  !<
  type :: LongLineReaderType

    character(len=:), allocatable :: line
    character(len=:), allocatable :: last_line
    integer(I4B) :: nbackspace = 0
    integer(I4B) :: iostat = 0
    integer(I4B) :: last_unit = 0

  contains

    procedure :: bkspc
    procedure :: rdcom

  end type LongLineReaderType

contains

  !> @brief Return the first non-comment line
  !!
  !! Skip through any comments and return the first
  !! non-commented line.  If an end of file was
  !! encountered previously, then return a blank line.
  !! If a backspace was called prior to this call,
  !! then do not read a new line and return last_line
  !! instead.
  !!
  !<
  subroutine rdcom(this, iu, iout, line, ierr)
    class(LongLineReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    character(len=:), intent(inout), allocatable :: line
    integer(I4B), intent(inout) :: ierr

    ierr = 0

    ! If using this reader to read from a new file
    ! then reset state
    if (iu /= this%last_unit) then
      this%nbackspace = 0
      this%iostat = 0
    end if

    if (this%nbackspace == 1) then
      ! If backspace was called, then return last line
      if (allocated(line)) deallocate (line)
      allocate (character(len=len(this%last_line) + 1) :: line)
      line(:) = this%last_line(:)
      this%nbackspace = 0
    else
      ! if end of file was reached previously, then return a
      ! blank line and return ierr as IOSTAT_END
      if (this%iostat == IOSTAT_END) then
        line = ' '
        ierr = IOSTAT_END
      else
        call u9rdcom(iu, iout, line, ierr)
      end if
      this%last_line = line
      this%iostat = ierr
    end if
    this%last_unit = iu
  end subroutine rdcom

  !> @brief Emulate a Fortran backspace
  !!
  !!  Emulate a fortran backspace call by storing
  !!  the current line in long_line
  !!
  !<
  subroutine bkspc(this, iin)
    class(LongLineReaderType) :: this
    integer(I4B), intent(in) :: iin
    if (this%nbackspace > 0) then
      call store_error( &
        "Programming error in LongLineReaderType%bkspc(). Backspace &
        & called more than once for an open file.", &
        terminate=.true.)
    else
      this%nbackspace = 1
    end if
  end subroutine bkspc

end module LongLineReaderModule

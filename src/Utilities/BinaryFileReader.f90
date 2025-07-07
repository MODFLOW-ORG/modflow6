module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream

  public :: BinaryFileReaderType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    integer(I4B) :: kstp
    integer(I4B) :: kper
    integer(I4B) :: kstpnext
    integer(I4B) :: kpernext
    logical(LGP) :: endoffile
    real(DP) :: delt
    real(DP) :: pertim
    real(DP) :: totim
  contains
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_if(this, success, iout)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
    end subroutine read_record_if
  end interface
contains

  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat

    if (.not. this%endoffile) then
      read (this%inunit, iostat=iostat) this%kstpnext, this%kpernext
      if (iostat == 0) then
        call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
      else if (iostat < 0) then
        this%endoffile = .true.
      end if
    end if
  end subroutine peek_record

end module BinaryFileReaderModule

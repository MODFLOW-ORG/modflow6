module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop

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
    procedure(read_record_inf), deferred :: read_record
    procedure(peek_record_inf), deferred :: peek_record
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_inf(this, success, iout, header_only)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
      logical(LGP), intent(in), optional :: header_only
    end subroutine read_record_inf
  end interface

  abstract interface
    subroutine peek_record_inf(this)
      import BinaryFileReaderType
      class(BinaryFileReaderType), intent(inout) :: this
    end subroutine peek_record_inf
  end interface
contains

end module BinaryFileReaderModule

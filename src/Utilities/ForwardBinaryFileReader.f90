module ForwardBinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use BinaryFileReaderModule, only: BinaryFileReaderType, &
                                    read_record_inf
  use InputOutputModule, only: fseek_stream

  private
  public :: ForwardBinaryFileReaderType

  type, extends(BinaryFileReaderType) :: ForwardBinaryFileReaderType
  contains
    procedure :: read_record_internal
    procedure :: read_record
    procedure :: peek_record
  end type ForwardBinaryFileReaderType

contains
  subroutine read_record(this, success, iout, header_only)
    class(ForwardBinaryFileReaderType), intent(inout) :: this
    logical(LGP), intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    logical(LGP), intent(in), optional :: header_only

    call this%read_record_internal(success, iout, header_only)
  end subroutine read_record

  subroutine read_record_internal(this, success, iout, header_only)
    class(ForwardBinaryFileReaderType), intent(inout) :: this
    logical(LGP), intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    logical(LGP), intent(in), optional :: header_only

    call pstop(1, 'Extending types must implement read_record_internal')
  end subroutine read_record_internal

  subroutine peek_record(this)
    class(ForwardBinaryFileReaderType), intent(inout) :: this
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
end module ForwardBinaryFileReaderModule

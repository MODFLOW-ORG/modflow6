module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream

  public :: BinaryFileHeaderType, &
            BinaryFileHeaderWrapperType, &
            BinaryFileReaderType

  type :: BinaryFileHeaderType
    integer(I4B) :: pos
    integer(I4B) :: kper, kstp
    real(DP) :: delt, pertim, totim
  contains
    procedure :: get_str
  end type BinaryFileHeaderType

  type :: BinaryFileHeaderWrapperType
    class(BinaryFileHeaderType), allocatable :: header
  end type BinaryFileHeaderWrapperType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    type(BinaryFileHeaderType) :: header
    type(BinaryFileHeaderType) :: headernext
    class(BinaryFileHeaderWrapperType), allocatable :: headers(:)
    integer(I4B) :: current
    integer(I4B) :: total
    logical(LGP) :: indexed
    logical(LGP) :: endoffile
  contains
    procedure(read_record_if), deferred :: read_record
    procedure :: peek_record
    procedure :: build_index
    procedure :: rewind
  end type BinaryFileReaderType

  abstract interface
    subroutine read_record_if(this, success, iout, header_only)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
      logical(LGP), intent(in), optional :: header_only
    end subroutine read_record_if
  end interface
contains

  !> @brief Get a string representation of the header.
  function get_str(this) result(str)
    class(BinaryFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str

    write (str, '(*(G0))') &
      'Binary file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ')'
    str = trim(str)
  end function get_str

  !> @brief Peek to see if another record is available.
  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat

    if (.not. this%endoffile) then
      read (this%inunit, iostat=iostat) this%headernext%kstp, this%headernext%kper
      if (iostat == 0) then
        call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
      else if (iostat < 0) then
        this%endoffile = .true.
      end if
    end if
  end subroutine peek_record

  subroutine build_index(this, iout)
    class(BinaryFileReaderType), intent(inout) :: this
    integer(I4B), intent(in), optional :: iout
    ! local
    integer(I4B) :: i
    logical(LGP) :: success

    if (this%indexed) return
    call this%rewind()
    this%total = 0
    i = 0
    do
      call this%read_record(success, iout, header_only=.true.)
      if (success) i = i + 1
      if (this%endoffile) exit
      if (.not. success) call pstop(1, 'Error reading record header')
    end do
    call this%rewind()
    this%total = i
    allocate (this%headers(this%total))
    i = 0
    do
      call this%read_record(success, iout, header_only=.true.)
      if (this%endoffile) exit
      if (.not. success) call pstop(1, 'Error reading record header')
      i = i + 1
      allocate (this%headers(i)%header, source=this%header)
    end do
    call this%rewind()
    this%indexed = .true.
  end subroutine build_index

  subroutine rewind(this)
    class(BinaryFileReaderType), intent(inout) :: this
    rewind(this%inunit)
    this%current = 0
    this%endoffile = .false.
  end subroutine rewind

end module BinaryFileReaderModule

module BinaryFileReaderModule

  use KindModule, only: I4B, I8B, DP, LGP
  use ErrorUtilModule, only: pstop
  use InputOutputModule, only: fseek_stream

  public :: BinaryFileHeaderType, BinaryFileReaderType

  type :: BinaryFileHeaderType
    integer(I4B) :: pos, size
    integer(I4B) :: kper, kstp
    real(DP) :: pertim, totim
  contains
    procedure :: get_str
  end type BinaryFileHeaderType

  type, abstract :: BinaryFileReaderType
    integer(I4B) :: inunit
    integer(I4B) :: nrecords
    logical(LGP) :: indexed
    logical(LGP) :: endoffile
    integer(I4B), allocatable :: record_sizes(:)
    class(BinaryFileHeaderType), allocatable :: header
    class(BinaryFileHeaderType), allocatable :: headernext
  contains
    procedure(read_header_if), deferred :: read_header
    procedure(read_record_if), deferred :: read_record
    procedure :: build_index
    procedure :: peek_record
    procedure :: rewind
  end type BinaryFileReaderType

  abstract interface
    subroutine read_header_if(this, success, iout)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
    end subroutine read_header_if

    subroutine read_record_if(this, success, iout)
      import BinaryFileReaderType
      import I4B, LGP
      class(BinaryFileReaderType), intent(inout) :: this
      logical(LGP), intent(out) :: success
      integer(I4B), intent(in), optional :: iout
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
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ')'
    str = trim(str)
  end function get_str

  !> @brief Build an index of data sizes in the file
  subroutine build_index(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    logical(LGP) :: success
    integer(I4B) :: pos, record_size

    this%indexed = .false.
    this%nrecords = 0
    call this%rewind()
    if (allocated(this%record_sizes)) deallocate (this%record_sizes)

    ! first pass: count
    do
      call this%read_record(success)
      if (.not. success) exit
      this%nrecords = this%nrecords + 1
    end do
    call this%rewind()
    allocate (this%record_sizes(this%nrecords))
    ! second pass: save
    do
      call this%read_record(success)
      inquire (this%inunit, pos=pos)
      record_size = pos - this%header%pos
      if (.not. success) exit
      this%record_sizes(this%nrecords) = this%header%size + record_size
    end do
    call this%rewind()

    this%indexed = .true.
  end subroutine build_index

  !> @brief Peek to see if another record is available.
  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat

    if (.not. this%endoffile) then
      inquire (this%inunit, pos=this%headernext%pos)
      read (this%inunit, iostat=iostat) this%headernext%kstp, this%headernext%kper
      if (iostat == 0) then
        call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
      else if (iostat < 0) then
        this%endoffile = .true.
        if (allocated(this%headernext)) deallocate (this%headernext)
      end if
    end if
  end subroutine peek_record

  !> @brief Rewind the file to the beginning.
  subroutine rewind (this)
    class(BinaryFileReaderType), intent(inout) :: this

    rewind (this%inunit)
    if (allocated(this%header)) deallocate (this%header)
    if (allocated(this%headernext)) deallocate (this%headernext)
    allocate (BinaryFileHeaderType :: this%header)
    allocate (BinaryFileHeaderType :: this%headernext)
    this%header%pos = 1
    this%headernext%pos = 1
    this%endoffile = .false.
  end subroutine rewind

end module BinaryFileReaderModule

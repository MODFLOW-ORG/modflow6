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
    logical(LGP) :: indexed = .false.
    logical(LGP) :: endoffile = .false.
    logical(LGP) :: backward = .false.
    integer(I4B) :: current_index = 0
    integer(I4B), allocatable :: record_positions(:)
    class(BinaryFileHeaderType), allocatable :: header
    class(BinaryFileHeaderType), allocatable :: headernext
  contains
    procedure(read_header_if), deferred :: read_header
    procedure(read_record_if), deferred :: read_record
    procedure :: build_index
    procedure :: peek_record
    procedure :: rewind
    procedure :: seek_to_index
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

  !> @brief Build an index of record positions in the file
  subroutine build_index(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    logical(LGP) :: success
    integer(I4B) :: i

    this%indexed = .false.
    this%nrecords = 0
    call this%rewind()
    if (allocated(this%record_positions)) deallocate (this%record_positions)

    ! first pass: count
    do
      call this%read_record(success)
      if (.not. success) exit
      this%nrecords = this%nrecords + 1
    end do

    call this%rewind()
    allocate (this%record_positions(this%nrecords))
    ! second pass: save positions (read first so header%pos is set by inquire)
    i = 0
    do
      call this%read_record(success)
      if (.not. success) exit
      i = i + 1
      this%record_positions(i) = this%header%pos
    end do
    call this%rewind()

    ! Initialize index position based on direction
    if (this%backward) then
      this%current_index = this%nrecords + 1
    else
      this%current_index = 0
    end if

    this%indexed = .true.
  end subroutine build_index

  !> @brief Peek to see if another record is available.
  subroutine peek_record(this)
    class(BinaryFileReaderType), intent(inout) :: this
    ! local
    integer(I4B) :: iostat, next_idx

    if (.not. this%endoffile) then
      if (this%indexed .and. this%backward) then
        ! Use index for backward peek
        if (this%backward) then
          next_idx = this%current_index - 1
        else
          next_idx = this%current_index + 1
        end if

        if (next_idx >= 1 .and. next_idx <= this%nrecords) then
          ! Seek to next position and read header
          call fseek_stream( &
            this%inunit, &
            this%record_positions(next_idx), &
            0, iostat)
          read (this%inunit, iostat=iostat) &
            this%headernext%kstp, this%headernext%kper
          ! Seek back to current position
          if (this%current_index >= 1 .and. &
              this%current_index <= this%nrecords) then
            call fseek_stream( &
              this%inunit, &
              this%record_positions(this%current_index), &
              0, iostat)
          end if
        else
          this%endoffile = .true.
        end if
      else
        ! Original forward peek logic
        inquire (this%inunit, pos=this%headernext%pos)
        read (this%inunit, iostat=iostat) &
          this%headernext%kstp, this%headernext%kper
        if (iostat == 0) then
          call fseek_stream(this%inunit, -2 * I4B, 1, iostat)
        else if (iostat < 0) then
          this%endoffile = .true.
          if (allocated(this%headernext)) deallocate (this%headernext)
        end if
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

    ! Reset index position based on direction
    if (this%indexed) then
      if (this%backward) then
        this%current_index = this%nrecords + 1
      else
        this%current_index = 0
      end if
    else
      this%current_index = 0
    end if
  end subroutine rewind

  !> @brief Seek to a specific index position
  subroutine seek_to_index(this, index)
    class(BinaryFileReaderType), intent(inout) :: this
    integer(I4B), intent(in) :: index
    ! local
    integer(I4B) :: iostat

    if (.not. this%indexed) then
      call pstop(1, 'Cannot seek to index in unindexed file')
    end if

    if (index < 1 .or. index > this%nrecords) then
      this%endoffile = .true.
      return
    end if

    call fseek_stream(this%inunit, this%record_positions(index), 0, iostat)
    if (iostat /= 0) then
      call pstop(1, 'Error seeking to index position')
    end if

    this%current_index = index
  end subroutine seek_to_index

end module BinaryFileReaderModule

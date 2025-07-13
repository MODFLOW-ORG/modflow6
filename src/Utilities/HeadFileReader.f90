module HeadFileReaderModule

  use KindModule
  use ConstantsModule, only: LINELENGTH
  use BinaryFileReaderModule, only: BinaryFileReaderType, BinaryFileHeaderType

  implicit none

  private
  public :: HeadFileReaderType, HeadFileHeaderType

  type, extends(BinaryFileHeaderType) :: HeadFileHeaderType
    character(len=16) :: text
    integer(I4B) :: ncol, nrow, ilay
  contains
    procedure :: get_str
  end type HeadFileHeaderType

  type, extends(BinaryFileReaderType) :: HeadFileReaderType
    integer(I4B) :: nlay
    real(DP), dimension(:), allocatable :: head
  contains
    procedure :: initialize
    procedure :: read_header
    procedure :: read_record
    procedure :: finalize
  end type HeadFileReaderType

contains

  !< @brief initialize
  !<
  subroutine initialize(this, iu, iout)
    ! -- dummy
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: kstp_last, kper_last
    logical :: success
    !
    this%inunit = iu
    this%endoffile = .false.
    this%nlay = 0
    !
    allocate (HeadFileHeaderType :: this%header)
    allocate (HeadFileHeaderType :: this%headernext)
    !
    ! -- Read the first head data record to set kstp_last, kstp_last
    call this%read_record(success)
    kstp_last = this%header%kstp
    kper_last = this%header%kper
    rewind (this%inunit)
    !
    ! -- Determine number of records within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading binary file to determine number of records per time step.'
    do
      call this%read_record(success, iout)
      if (.not. success) exit
      if (kstp_last /= this%header%kstp .or. kper_last /= this%header%kper) exit
      this%nlay = this%nlay + 1
    end do
    rewind (this%inunit)
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nlay, &
      ' unique records in binary file.'
  end subroutine initialize

  !< @brief read header only
  !<
  subroutine read_header(this, success, iout)
    ! -- dummy
    class(HeadFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! -- local
    integer(I4B) :: iostat
    !
    success = .true.
    select type (h => this%header)
    type is (HeadFileHeaderType)
      h%kstp = 0
      h%kper = 0
      h%text = ''
      h%ncol = 0
      h%nrow = 0
      h%ilay = 0
      read (this%inunit, iostat=iostat) h%kstp, h%kper, &
        h%pertim, h%totim, h%text, h%ncol, h%nrow, h%ilay
      if (iostat /= 0) then
        success = .false.
        if (iostat < 0) this%endoffile = .true.
        return
      end if
    end select
  end subroutine read_header

  !< @brief read record
  !<
  subroutine read_record(this, success, iout)
    ! -- modules
    use InputOutputModule, only: fseek_stream
    ! -- dummy
    class(HeadFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    ! -- local
    integer(I4B) :: iout_opt
    integer(I4B) :: ncol, nrow
    !
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    !
    call this%read_header(success, iout_opt)
    if (.not. success) return
    !
    select type (h => this%header)
    type is (HeadFileHeaderType)
      ncol = h%ncol
      nrow = h%nrow
    end select
    !
    ! -- allocate head to proper size
    if (.not. allocated(this%head)) then
      allocate (this%head(ncol * nrow))
    else
      if (size(this%head) /= ncol * nrow) then
        deallocate (this%head)
        allocate (this%head(ncol * nrow))
      end if
    end if
    !
    ! -- read the head array
    read (this%inunit) this%head
    !
    call this%peek_record()
  end subroutine read_record

  !< @brief finalize
  !<
  subroutine finalize(this)
    class(HeadFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%head)) deallocate (this%head)
    if (allocated(this%header)) deallocate (this%header)
    if (allocated(this%headernext)) deallocate (this%headernext)
  end subroutine finalize

  !> @brief Get a string representation of the head file header.
  function get_str(this) result(str)
    class(HeadFileHeaderType), intent(in) :: this
    character(len=:), allocatable :: str

    write (str, '(*(G0))') &
      'Head file header (pos: ', this%pos, &
      ', kper: ', this%kper, &
      ', kstp: ', this%kstp, &
      ', delt: ', this%delt, &
      ', pertim: ', this%pertim, &
      ', totim: ', this%totim, &
      ', text: ', trim(this%text), &
      ', ncol: ', this%ncol, &
      ', nrow: ', this%nrow, &
      ', ilay: ', this%ilay, &
      ')'
    str = trim(str)
  end function get_str

end module HeadFileReaderModule

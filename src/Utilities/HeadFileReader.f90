module HeadFileReaderModule

  use KindModule
  use ConstantsModule, only: LINELENGTH
  use BinaryFileReaderModule, only: BinaryFileReaderType
  use InputOutputModule, only: fseek_stream

  implicit none

  private
  public :: HeadFileReaderType

  type, extends(BinaryFileReaderType) :: HeadFileReaderType
    character(len=16) :: text
    integer(I4B) :: nlay
    real(DP), dimension(:), allocatable :: head
  contains
    procedure :: initialize
    procedure :: read_record
    procedure :: finalize
  end type HeadFileReaderType

contains

  !< @brief initialize
  !<
  subroutine initialize(this, iu, iout)
    ! dummy
    class(HeadFileReaderType) :: this
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    ! local
    integer(I4B) :: kstp_last, kper_last
    logical :: success
    
    this%inunit = iu
    this%endoffile = .false.
    this%nlay = 0
    
    call this%build_index(iout)
    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Indexed ', this%total, &
      ' total records in head file.'
    
    ! Read the first head data record to set kstp_last, kstp_last
    call this%read_record(success, header_only=.true.)
    kstp_last = this%header%kstp
    kper_last = this%header%kper
    call this%rewind()
    
    ! Determine number of records within a time step
    if (iout > 0) &
      write (iout, '(a)') &
      'Reading binary file to determine number of records per time step.'
    do
      call this%read_record(success, iout, header_only=.true.)
      if (.not. success) exit
      if (kstp_last /= this%header%kstp .or. kper_last /= this%header%kper) exit
      this%nlay = this%nlay + 1
    end do
    call this%rewind()

    if (iout > 0) &
      write (iout, '(a, i0, a)') 'Detected ', this%nlay, &
      ' unique records in binary file.'
  end subroutine initialize

  !< @brief read record
  !<
  subroutine read_record(this, success, iout, header_only)
    ! dummy
    class(HeadFileReaderType), intent(inout) :: this
    logical, intent(out) :: success
    integer(I4B), intent(in), optional :: iout
    logical, intent(in), optional :: header_only
    ! local
    logical(LGP) :: header_only_opt
    integer(I4B) :: iostat, iout_opt
    integer(I4B) :: ncol, nrow, ilay
    
    if (present(header_only)) then
      header_only_opt = header_only
    else
      header_only_opt = .false.
    end if
    if (present(iout)) then
      iout_opt = iout
    else
      iout_opt = 0
    end if
    
    this%header%kstp = 0
    this%header%kper = 0
    success = .true.
    this%headernext%kstp = 0
    this%headernext%kper = 0
    inquire (unit=this%inunit, pos=this%header%pos)
    read (this%inunit, iostat=iostat) this%header%kstp, this%header%kper, &
      this%header%pertim, this%header%totim, this%text, ncol, nrow, ilay
    if (iostat /= 0) then
      success = .false.
      if (iostat < 0) this%endoffile = .true.
      return
    end if
    this%current = this%current + 1
    if (header_only_opt) then
      call fseek_stream(this%inunit, ncol * nrow * 8, 1, iostat)
      call this%peek_record()
      return
    end if
    
    ! allocate head to proper size
    if (.not. allocated(this%head)) then
      allocate (this%head(ncol * nrow))
    else
      if (size(this%head) /= ncol * nrow) then
        deallocate (this%head)
        allocate (this%head(ncol * nrow))
      end if
    end if
    
    ! read the head array
    read (this%inunit) this%head
    
    call this%peek_record()
  end subroutine read_record

  !< @brief finalize
  !<
  subroutine finalize(this)
    class(HeadFileReaderType) :: this
    close (this%inunit)
    if (allocated(this%head)) deallocate (this%head)
  end subroutine finalize

end module HeadFileReaderModule

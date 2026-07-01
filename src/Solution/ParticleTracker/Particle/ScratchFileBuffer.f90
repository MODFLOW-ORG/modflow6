!> @brief Scratch-file particle event buffer module.
!<
module ScratchFileBufferModule

  use KindModule, only: I4B
  use ErrorUtilModule, only: pstop
  use ParticleTrackEventBufferModule, only: ParticleTrackEventBufferType, &
                                            ParticleTrackRecordType, &
                                            ParticleTrackFileType, &
                                            save_record

  implicit none
  private
  public :: ScratchFileBufferType

  !> @brief Scratch-file particle event buffer. Records are written to
  !! an unformatted sequential scratch file that is rewound on discard
  !! and read back sequentially on flush.
  type, extends(ParticleTrackEventBufferType) :: ScratchFileBufferType
    integer(I4B) :: iun = 0 !< scratch file unit
  contains
    procedure :: init => scratch_init
    procedure :: append => scratch_append
    procedure :: flush => scratch_flush
    procedure :: discard => scratch_discard
    procedure :: destroy => scratch_destroy
  end type ScratchFileBufferType

contains

  subroutine scratch_init(this)
    class(ScratchFileBufferType) :: this
    integer(I4B) :: istat
    open (newunit=this%iun, status='scratch', form='unformatted', &
          access='sequential', iostat=istat)
    if (istat /= 0) &
      call pstop(1, 'failed to open scratch track buffer file')
  end subroutine scratch_init

  subroutine scratch_append(this, rec)
    class(ScratchFileBufferType) :: this
    type(ParticleTrackRecordType), intent(in) :: rec
    write (this%iun) rec
    this%nrecords = this%nrecords + 1
  end subroutine scratch_append

  subroutine scratch_flush(this, files)
    class(ScratchFileBufferType) :: this
    type(ParticleTrackFileType), intent(in) :: files(:)
    integer(I4B) :: n, i
    type(ParticleTrackRecordType) :: rec

    rewind (this%iun)
    do n = 1, this%nrecords
      read (this%iun) rec
      do i = 1, size(files)
        if (files(i)%iun > 0 .and. &
            (files(i)%iprp == -1 .or. files(i)%iprp == rec%iprp)) &
          call save_record(files(i)%iun, rec, files(i)%csv)
      end do
    end do
    rewind (this%iun)
    this%nrecords = 0
  end subroutine scratch_flush

  subroutine scratch_discard(this)
    class(ScratchFileBufferType) :: this
    rewind (this%iun)
    this%nrecords = 0
  end subroutine scratch_discard

  subroutine scratch_destroy(this)
    class(ScratchFileBufferType) :: this
    if (this%iun > 0) then
      close (this%iun)
      this%iun = 0
    end if
    this%nrecords = 0
  end subroutine scratch_destroy

end module ScratchFileBufferModule

!> @brief In-memory particle event buffer module.
!<
module MemoryBufferModule

  use KindModule, only: I4B
  use ParticleTrackEventBufferModule, only: ParticleTrackEventBufferType, &
                                            TrackRecordType, &
                                            ParticleTrackFileType, &
                                            save_record

  implicit none
  private
  public :: MemoryBufferType

  !> @brief In-memory particle event buffer. Records are held in a
  !! dynamically growing array that doubles in capacity as needed.
  type, extends(ParticleTrackEventBufferType) :: MemoryBufferType
    type(TrackRecordType), allocatable :: records(:) !< buffer
  contains
    procedure :: append => memory_append
    procedure :: flush => memory_flush
    procedure :: discard => memory_discard
    procedure :: destroy => memory_destroy
  end type MemoryBufferType

contains

  subroutine memory_append(this, rec)
    class(MemoryBufferType) :: this
    type(TrackRecordType), intent(in) :: rec
    type(TrackRecordType), allocatable :: tmp(:)

    if (.not. allocated(this%records)) then
      allocate (this%records(64))
    else if (this%nrecords == size(this%records)) then
      allocate (tmp(size(this%records) * 2))
      tmp(1:this%nrecords) = this%records
      call move_alloc(tmp, this%records)
    end if
    this%nrecords = this%nrecords + 1
    this%records(this%nrecords) = rec
  end subroutine memory_append

  subroutine memory_flush(this, files)
    class(MemoryBufferType) :: this
    type(ParticleTrackFileType), intent(in) :: files(:)
    integer(I4B) :: n, i
    type(TrackRecordType) :: rec

    do n = 1, this%nrecords
      rec = this%records(n)
      do i = 1, size(files)
        if (files(i)%iun > 0 .and. &
            (files(i)%iprp == -1 .or. files(i)%iprp == rec%iprp)) &
          call save_record(files(i)%iun, rec, files(i)%csv)
      end do
    end do
    this%nrecords = 0
  end subroutine memory_flush

  subroutine memory_discard(this)
    class(MemoryBufferType) :: this
    this%nrecords = 0
  end subroutine memory_discard

  subroutine memory_destroy(this)
    class(MemoryBufferType) :: this
    if (allocated(this%records)) deallocate (this%records)
  end subroutine memory_destroy

end module MemoryBufferModule

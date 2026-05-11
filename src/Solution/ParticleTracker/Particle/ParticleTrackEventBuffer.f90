!> @brief Particle track event buffer module.
!!
!! Defines the abstract ParticleTrackEventBufferType, the TrackRecordType
!! and ParticleTrackFileType data types, and the save_record helper used
!! by concrete buffer implementations.
!<
module ParticleTrackEventBufferModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop

  implicit none

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

  !> @brief Flat record of a particle track event.
  type :: TrackRecordType
    integer(I4B) :: kper, kstp, imdl, iprp, irpt, ilay, icu, izone
    integer(I4B) :: istatus, ireason
    real(DP) :: trelease, ttrack, x, y, z
    character(len=40) :: name
  end type TrackRecordType

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! Can be associated with a particle release point (PRP) package
  !! or with an entire model, and can be binary or text (CSV).
  !<
  type :: ParticleTrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type ParticleTrackFileType

  !> @brief Event buffering strategy
  type, abstract :: ParticleTrackEventBufferType
    integer(I4B) :: nrecords = 0 !< number of records stored
  contains
    procedure :: init => buffer_init !< open/allocate resources
    procedure(buffer_append), deferred :: append !< buffer one record
    procedure(buffer_flush), deferred :: flush !< write buffered records, reset
    procedure(buffer_simple), deferred :: discard !< reset without writing
    procedure(buffer_simple), deferred :: destroy !< release resources
  end type ParticleTrackEventBufferType

  abstract interface
    subroutine buffer_append(this, rec)
      import ParticleTrackEventBufferType, TrackRecordType
      class(ParticleTrackEventBufferType) :: this
      type(TrackRecordType), intent(in) :: rec
    end subroutine buffer_append

    subroutine buffer_flush(this, files)
      import ParticleTrackEventBufferType, ParticleTrackFileType
      class(ParticleTrackEventBufferType) :: this
      type(ParticleTrackFileType), intent(in) :: files(:)
    end subroutine buffer_flush

    ! Shared interface for discard and destroy: both take only `this`.
    subroutine buffer_simple(this)
      import ParticleTrackEventBufferType
      class(ParticleTrackEventBufferType) :: this
    end subroutine buffer_simple
  end interface

contains

  subroutine buffer_init(this)
    class(ParticleTrackEventBufferType) :: this
  end subroutine buffer_init

  !> @brief Save an event record to a binary or CSV file.
  subroutine save_record(iun, rec, csv)
    integer(I4B), intent(in) :: iun
    type(TrackRecordType), intent(in) :: rec
    logical(LGP), intent(in) :: csv

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, trim(rec%name)
    else
      write (iun) &
        rec%kper, rec%kstp, rec%imdl, rec%iprp, rec%irpt, &
        rec%ilay, rec%icu, rec%izone, rec%istatus, rec%ireason, &
        rec%trelease, rec%ttrack, rec%x, rec%y, rec%z, rec%name
    end if
  end subroutine save_record

end module ParticleTrackEventBufferModule

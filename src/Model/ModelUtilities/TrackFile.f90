module TrackFileModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DPIO180
  use ParticleModule, only: ParticleType, ACTIVE
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform

  implicit none
  public :: TrackFileType
  public :: save_record

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! Can be associated with a particle release point (PRP) package
  !! or with an entire model, and can be binary or comma-separated.
  !!
  !! Each particle's pathline consists of 1+ records reported as the particle
  !! is tracked over the model domain. Records are snapshots of the particle's
  !! state (e.g. tracking status, position) at a particular moment in time.
  !!
  !! Particles have no ID property. Particles can be uniquely identified
  !! by composite key, i.e. combination of fields:
  !!
  !!   - imdl: originating model ID
  !!   - iprp: originating PRP ID
  !!   - irpt: particle release location ID
  !!   - trelease: particle release time
  !<
  type :: TrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type TrackFileType

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

contains

  !> @brief Save a particle track record to a binary or CSV file.
  subroutine save_record(iun, particle, dis, kper, kstp, reason, csv)
    ! dummy
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: reason
    logical(LGP), intent(in) :: csv
    ! local
    real(DP) :: x, y, z, xout, yout, zout
    integer(I4B) :: status

    ! Convert from cell-local to model coordinates if needed
    call particle%get_model_coords(x, y, z)

    ! Apply model grid offset/rotation if needed
    if (particle%icoords == 1 .and. &
        (dis%angrot /= DZERO .or. &
         dis%xorigin /= DZERO .or. &
         dis%yorigin /= DZERO)) then
      call transform(x, y, z, &
                     xout, yout, zout, &
                     dis%xorigin, &
                     dis%yorigin, &
                     DZERO, &
                     sin(dis%angrot * DPIO180), &
                     cos(dis%angrot * DPIO180), &
                     invert=.true.)
      x = xout
      y = yout
    end if

    ! Set status
    if (particle%istatus .lt. 0) then
      status = ACTIVE
    else
      status = particle%istatus
    end if

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        kper, &
        kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        reason, &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        trim(adjustl(particle%name))
    else
      write (iun) &
        kper, &
        kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        reason, &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        particle%name
    end if
  end subroutine
end module TrackFileModule

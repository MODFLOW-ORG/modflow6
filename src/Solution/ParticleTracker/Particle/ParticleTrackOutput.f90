module ParticleTrackOutputModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, DPIO180
  use ParticleModule, only: ParticleType, ACTIVE
  use ParticleEventModule, only: ParticleEventType
  use ParticleEventsModule, only: ParticleEventConsumerType, &
                                  ParticleEventDispatcherType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform

  implicit none
  public :: ParticleTrackFileType, &
            ParticleTrackOutputType, &
            ParticleTrackEventSelectionType
  private :: log_event, save_event

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
  type :: ParticleTrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type ParticleTrackFileType

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

  type :: ParticleTrackEventSelectionType
    logical(LGP) :: release !< track release events
    logical(LGP) :: cellexit !< track cell exits
    logical(LGP) :: timestep !< track timestep ends
    logical(LGP) :: terminate !< track termination events
    logical(LGP) :: weaksink !< track weak sink exit events
    logical(LGP) :: usertime !< track user-selected times
  end type ParticleTrackEventSelectionType

  !> @brief Manages particle track (i.e. pathline) output (logging/writing).
  !!
  !! Optionally filters events as selected in the PRT Output Control package.
  !! An arbitrary number of files can be managed, resizing is done as needed.
  !<
  type, extends(ParticleEventConsumerType) :: ParticleTrackOutputType
    private
    integer(I4B), public :: iout = -1 !<  log file unit
    integer(I4B), public :: ntrackfiles !< number of track files
    type(ParticleTrackFileType), public, allocatable :: trackfiles(:) !< track files
    type(ParticleTrackEventSelectionType), public :: selection !< event selection
  contains
    procedure, public :: init_track_file
    procedure, public :: destroy
    procedure, public :: select
    procedure :: expand
    procedure :: handle_event
    procedure :: is_selected
    procedure :: should_save
    procedure :: should_log
  end type ParticleTrackOutputType

contains

  !> @brief Initialize a new track file
  subroutine init_track_file(this, iun, csv, iprp)
    ! dummy
    class(ParticleTrackOutputType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    ! local
    type(ParticleTrackFileType), pointer :: file

    ! Allocate or expand array
    if (.not. allocated(this%trackfiles)) then
      allocate (this%trackfiles(1))
    else
      call this%expand(increment=1)
    end if

    ! Setup new file
    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp

    ! Update array and counter
    this%ntrackfiles = size(this%trackfiles)
    this%trackfiles(this%ntrackfiles) = file

  end subroutine init_track_file

  subroutine destroy(this)
    class(ParticleTrackOutputType) :: this

    if (allocated(this%trackfiles)) deallocate (this%trackfiles)
  end subroutine destroy

  !> @brief Expand the trackfile array
  subroutine expand(this, increment)
    ! dummy
    class(ParticleTrackOutputType) :: this
    integer(I4B), optional, intent(in) :: increment
    ! local
    integer(I4B) :: inclocal
    integer(I4B) :: isize
    integer(I4B) :: newsize
    type(ParticleTrackFileType), allocatable, dimension(:) :: temp

    ! Initialize optional args
    if (present(increment)) then
      inclocal = increment
    else
      inclocal = 1
    end if

    ! Increase size of array
    if (allocated(this%trackfiles)) then
      isize = size(this%trackfiles)
      newsize = isize + inclocal
      allocate (temp(newsize))
      temp(1:isize) = this%trackfiles
      deallocate (this%trackfiles)
      call move_alloc(temp, this%trackfiles)
    else
      allocate (this%trackfiles(inclocal))
    end if

  end subroutine expand

  !> @brief Select events.
  subroutine select(this, &
                    release, &
                    cellexit, &
                    timestep, &
                    terminate, &
                    weaksink, &
                    usertime)
    class(ParticleTrackOutputType) :: this
    logical(LGP), intent(in) :: release
    logical(LGP), intent(in) :: cellexit
    logical(LGP), intent(in) :: timestep
    logical(LGP), intent(in) :: terminate
    logical(LGP), intent(in) :: weaksink
    logical(LGP), intent(in) :: usertime
    this%selection%release = release
    this%selection%cellexit = cellexit
    this%selection%timestep = timestep
    this%selection%terminate = terminate
    this%selection%weaksink = weaksink
    this%selection%usertime = usertime
  end subroutine select

  !> @brief Check if the given event code is selected.
  logical function is_selected(this, event_code) result(selected)
    class(ParticleTrackOutputType), intent(inout) :: this
    integer(I4B), intent(in) :: event_code

    selected = (this%selection%release .and. event_code == 0) .or. &
               (this%selection%cellexit .and. event_code == 1) .or. &
               (this%selection%timestep .and. event_code == 2) .or. &
               (this%selection%terminate .and. event_code == 3) .or. &
               (this%selection%weaksink .and. event_code == 4) .or. &
               (this%selection%usertime .and. event_code == 5)
  end function is_selected

  !> @brief Check whether a particle belongs in a given file, i.e.
  !! if the file is enabled and its group matches the particle's.
  logical function should_save(this, particle, file) result(save)
    class(ParticleTrackOutputType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    type(ParticleTrackFileType), intent(in) :: file

    save = (file%iun > 0 .and. &
            (file%iprp == -1 .or. file%iprp == particle%iprp))
  end function should_save

  !> @brief Save an event to a binary or CSV file.
  subroutine save_event(iun, particle, event, csv)
    ! dummy
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP), intent(in) :: csv
    ! local
    real(DP) :: x, y, z
    integer(I4B) :: status

    ! Convert from cell-local to model coordinates if needed
    call particle%get_model_coords(x, y, z)

    ! Set status
    if (particle%istatus .lt. 0) then
      status = ACTIVE
    else
      status = particle%istatus
    end if

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        event%kper, &
        event%kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        event%get_code(), &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        trim(adjustl(particle%name))
    else
      write (iun) &
        event%kper, &
        event%kstp, &
        particle%imdl, &
        particle%iprp, &
        particle%irpt, &
        particle%ilay, &
        particle%icu, &
        particle%izone, &
        status, &
        event%get_code(), &
        particle%trelease, &
        particle%ttrack, &
        x, &
        y, &
        z, &
        particle%name
    end if
  end subroutine save_event

  logical function should_log(this)
    class(ParticleTrackOutputType), intent(inout) :: this

    should_log = this%iout >= 0
  end function should_log

  !> @brief Print a particle event summary.
  subroutine log_event(iun, particle, event)
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event

    if (iun >= 0) &
      write (iun, '(*(G0))') &
      'Particle (Model: ', particle%imdl, &
      '  PRP Package: ', particle%iprp, &
      '  Release Point: ', particle%irpt, &
      '  Release Time: ', particle%trelease, &
      '  Name: ', trim(adjustl(particle%name)), &
      ') ', event%get_str(), &
      ' in (Layer: ', particle%ilay, &
      '  Cell: ', particle%icu, &
      '  Zone: ', particle%izone, &
      ') at (Time ', particle%ttrack, &
      '  Period ', event%kper, &
      '  Timestep ', event%kstp, &
      ') with (Status: ', particle%istatus, &
      '  Context: ', trim(adjustl(event%context)), ')'
  end subroutine log_event

  !> @brief Handle a particle event.
  subroutine handle_event(this, particle, event)
    ! dummy
    class(ParticleTrackOutputType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    integer(I4B) :: i
    type(ParticleTrackFileType) :: file

    if (this%should_log()) &
      call log_event(this%iout, particle, event)

    if (this%is_selected(event%get_code())) then
      do i = 1, this%ntrackfiles
        file = this%trackfiles(i)
        if (this%should_save(particle, file)) &
          call save_event(file%iun, particle, event, csv=file%csv)
      end do
    end if
  end subroutine handle_event

end module ParticleTrackOutputModule

!> @brief Particle track output module.
!!
!! Each particle's track consists of events reported as the particle is
!! advected through the model domain. Events are snapshots of particle
!! state, along with optional metadata, at a particular moment in time.
!!
!! Particles have no ID property. A particle can be uniquely identified
!! by the unique combination of its release attributes (model, package,
!! position, and time). This is possible because only one particle may
!! be released from a given point at a given time.
!!
!! This module consumes particle events and is responsible for writing
!! them to one or more track files, binary or CSV, and for logging the
!! events if requested. Each track file is associated with either a PRP
!! package or with the full PRT model (there may only be 1 such latter).
!<
module ParticleTracksModule

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
            ParticleTracksType, &
            ParticleTrackEventSelectionType
  private :: log_event, save_event

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

  !> @brief Output file containing all or some particle pathlines.
  !!
  !! Can be associated with a particle release point (PRP) package
  !! or with an entire model, and can be binary or comma-separated.
  !<
  type :: ParticleTrackFileType
    private
    integer(I4B), public :: iun = 0 !< file unit number
    logical(LGP), public :: csv = .false. !< whether the file is binary or CSV
    integer(I4B), public :: iprp = -1 !< -1 is model-level file, 0 is exchange PRP
  end type ParticleTrackFileType

  !> @brief Selection of particle events.
  type :: ParticleTrackEventSelectionType
    logical(LGP) :: release !< track release events
    logical(LGP) :: cellexit !< track cell exits
    logical(LGP) :: timestep !< track timestep ends
    logical(LGP) :: terminate !< track termination events
    logical(LGP) :: weaksink !< track weak sink exit events
    logical(LGP) :: usertime !< track user-selected times
  end type ParticleTrackEventSelectionType

  !> @brief Manages particle track output (logging/writing).
  !!
  !! Optionally filters events as selected in the PRT Output Control package.
  !! An arbitrary number of files can be managed, resizing is done as needed.
  !<
  type, extends(ParticleEventConsumerType) :: ParticleTracksType
    private
    integer(I4B), public :: iout = -1 !<  log file unit
    integer(I4B), public :: ntrackfiles !< number of track files
    type(ParticleTrackFileType), public, allocatable :: files(:) !< track files
    type(ParticleTrackEventSelectionType), public :: selected !< event selection
  contains
    procedure, public :: init_file
    procedure, public :: is_selected
    procedure, public :: select_events
    procedure, public :: destroy
    procedure :: expand_files
    procedure :: handle_event
    procedure :: should_save
    procedure :: should_log
  end type ParticleTracksType

contains

  !> @brief Initialize a binary or CSV file.
  subroutine init_file(this, iun, csv, iprp)
    ! dummy
    class(ParticleTracksType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    ! local
    type(ParticleTrackFileType), pointer :: file

    if (.not. allocated(this%files)) then
      allocate (this%files(1))
    else
      call this%expand_files(increment=1)
    end if

    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp
    this%ntrackfiles = size(this%files)
    this%files(this%ntrackfiles) = file
  end subroutine init_file

  subroutine destroy(this)
    class(ParticleTracksType) :: this
    if (allocated(this%files)) deallocate (this%files)
  end subroutine destroy

  !> @brief Grow the array of track files.
  subroutine expand_files(this, increment)
    ! dummy
    class(ParticleTracksType) :: this
    integer(I4B), optional, intent(in) :: increment
    ! local
    integer(I4B) :: inclocal
    integer(I4B) :: isize
    integer(I4B) :: newsize
    type(ParticleTrackFileType), allocatable, dimension(:) :: temp

    if (present(increment)) then
      inclocal = increment
    else
      inclocal = 1
    end if

    if (allocated(this%files)) then
      isize = size(this%files)
      newsize = isize + inclocal
      allocate (temp(newsize))
      temp(1:isize) = this%files
      deallocate (this%files)
      call move_alloc(temp, this%files)
    else
      allocate (this%files(inclocal))
    end if
  end subroutine expand_files

  !> @brief Pick events to track.
  subroutine select_events(this, &
                           release, &
                           cellexit, &
                           timestep, &
                           terminate, &
                           weaksink, &
                           usertime)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: release
    logical(LGP), intent(in) :: cellexit
    logical(LGP), intent(in) :: timestep
    logical(LGP), intent(in) :: terminate
    logical(LGP), intent(in) :: weaksink
    logical(LGP), intent(in) :: usertime
    this%selected%release = release
    this%selected%cellexit = cellexit
    this%selected%timestep = timestep
    this%selected%terminate = terminate
    this%selected%weaksink = weaksink
    this%selected%usertime = usertime
  end subroutine select_events

  !> @brief Check if a given event code is selected for tracking.
  logical function is_selected(this, event_code) result(selected)
    class(ParticleTracksType), intent(inout) :: this
    integer(I4B), intent(in) :: event_code

    selected = (this%selected%release .and. event_code == 0) .or. &
               (this%selected%cellexit .and. event_code == 1) .or. &
               (this%selected%timestep .and. event_code == 2) .or. &
               (this%selected%terminate .and. event_code == 3) .or. &
               (this%selected%weaksink .and. event_code == 4) .or. &
               (this%selected%usertime .and. event_code == 5)
  end function is_selected

  !> @brief Check whether a particle belongs in a given file i.e.
  !! if the file is enabled and its group matches the particle's.
  logical function should_save(this, particle, file) result(save)
    class(ParticleTracksType), intent(inout) :: this
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

  !> @brief Log output unit valid?
  logical function should_log(this)
    class(ParticleTracksType), intent(inout) :: this
    should_log = this%iout >= 0
  end function should_log

  !> @brief Print a particle event summary.
  subroutine log_event(iun, particle, event)
    integer(I4B), intent(in) :: iun
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    character(len=:), allocatable :: particlename

    particlename = trim(adjustl(particle%name))
    if (len_trim(particlename) == 0) particlename = 'anonymous'

    if (iun >= 0) &
      write (iun, '(*(G0))') &
      'Particle (Model: ', particle%imdl, &
      ', Package: ', particle%iprp, &
      ', Point: ', particle%irpt, ' [', particlename, ']', &
      ', Time: ', particle%trelease, &
      ') ', event%get_str(), &
      ' in (Layer: ', particle%ilay, &
      ', Cell: ', particle%icu, &
      ', Zone: ', particle%izone, &
      ') at (X: ', particle%x, &
      ', Y: ', particle%y, &
      ', Z: ', particle%z, &
      ', Time: ', particle%ttrack, &
      ', Period: ', event%kper, &
      ', Timestep: ', event%kstp, &
      ') with (Status: ', particle%istatus, ')'
  end subroutine log_event

  !> @brief Handle a particle event.
  subroutine handle_event(this, particle, event)
    ! dummy
    class(ParticleTracksType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    ! local
    integer(I4B) :: i
    type(ParticleTrackFileType) :: file

    if (this%should_log()) &
      call log_event(this%iout, particle, event)

    if (this%is_selected(event%get_code())) then
      do i = 1, this%ntrackfiles
        file = this%files(i)
        if (this%should_save(particle, file)) &
          call save_event(file%iun, particle, event, csv=file%csv)
      end do
    end if
  end subroutine handle_event

end module ParticleTracksModule

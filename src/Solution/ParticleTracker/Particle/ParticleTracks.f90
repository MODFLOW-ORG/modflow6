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
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, DONE, DPIO180
  use ParticleModule, only: ParticleType, ACTIVE
  use ParticleEventModule, only: ParticleEventType
  use ReleaseEventModule, only: ReleaseEventType
  use TimeStepEventModule, only: TimeStepEventType
  use TerminationEventModule, only: TerminationEventType
  use WeakSinkEventModule, only: WeakSinkEventType
  use UserTimeEventModule, only: UserTimeEventType
  use ModelExitEventModule, only: ModelExitEventType
  use CellExitEventModule, only: CellExitEventType
  use SubcellExitEventModule, only: SubcellExitEventType
  use DroppedEventModule, only: DroppedEventType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform

  implicit none
  public :: ParticleTrackFileType, &
            ParticleTracksType, &
            ParticleTrackEventSelectionType, &
            write_particle_event
  private :: save_event

  character(len=*), parameter, public :: TRACKHEADER = &
    'kper,kstp,imdl,iprp,irpt,ilay,icell,izone,&
    &istatus,ireason,trelease,t,x,y,z,name'

  character(len=*), parameter, public :: TRACKDTYPES = &
    '<i4,<i4,<i4,<i4,<i4,<i4,<i4,<i4,&
    &<i4,<i4,<f8,<f8,<f8,<f8,<f8,|S40'

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

  !> @brief Selection of particle events.
  type :: ParticleTrackEventSelectionType
    logical(LGP) :: release !< track release events
    logical(LGP) :: featexit !< track grid feature exits
    logical(LGP) :: timestep !< track timestep ends
    logical(LGP) :: terminate !< track termination events
    logical(LGP) :: weaksink !< track weak sink exit events
    logical(LGP) :: usertime !< track user-selected times
    logical(LGP) :: subfexit !< track subfeature exits
    logical(LGP) :: dropped !< track water table drops
  end type ParticleTrackEventSelectionType

  !> @brief Particle track output manager. Handles printing as well as writing
  !! to files. One output unit can be configured for printing. Multiple files
  !! can be configured for writing, with each file optionally associated with
  !! a PRP package or with the full model. Events can be filtered by type, so
  !! that only certain event types are printed or written to files.
  !<
  type :: ParticleTracksType
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
    procedure :: should_save
    procedure :: should_print
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

  !> @brief Destroy the particle track manager.
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
                           featexit, &
                           timestep, &
                           terminate, &
                           weaksink, &
                           usertime, &
                           subfexit, &
                           dropped)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: release
    logical(LGP), intent(in) :: featexit
    logical(LGP), intent(in) :: timestep
    logical(LGP), intent(in) :: terminate
    logical(LGP), intent(in) :: weaksink
    logical(LGP), intent(in) :: usertime
    logical(LGP), intent(in) :: subfexit
    logical(LGP), intent(in) :: dropped
    this%selected%release = release
    this%selected%featexit = featexit
    this%selected%timestep = timestep
    this%selected%terminate = terminate
    this%selected%weaksink = weaksink
    this%selected%usertime = usertime
    this%selected%subfexit = subfexit
    this%selected%dropped = dropped
  end subroutine select_events

  !> @brief Check if a given event code is selected for tracking.
  logical function is_selected(this, event) result(selected)
    class(ParticleTracksType), intent(inout) :: this
    class(ParticleEventType), intent(in) :: event

    select type (event)
    type is (ReleaseEventType)
      selected = this%selected%release
    type is (ModelExitEventType)
      selected = this%selected%featexit
    type is (CellExitEventType)
      selected = this%selected%featexit
    type is (SubcellExitEventType)
      selected = this%selected%subfexit
    type is (TimeStepEventType)
      selected = this%selected%timestep
    type is (TerminationEventType)
      selected = this%selected%terminate
    type is (WeakSinkEventType)
      selected = this%selected%weaksink
    type is (UserTimeEventType)
      selected = this%selected%usertime
    type is (DroppedEventType)
      selected = this%selected%dropped
    class default
      call pstop(1, "unknown event type")
      selected = .false.
    end select

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

    if (csv) then
      write (iun, '(*(G0,:,","))') &
        event%kper, &
        event%kstp, &
        event%imdl, &
        event%iprp, &
        event%irpt, &
        event%ilay, &
        event%icu, &
        event%izone, &
        event%istatus, &
        event%get_code(), &
        event%trelease, &
        event%ttrack, &
        event%x, &
        event%y, &
        event%z, &
        trim(adjustl(particle%name))
    else
      write (iun) &
        event%kper, &
        event%kstp, &
        event%imdl, &
        event%iprp, &
        event%irpt, &
        event%ilay, &
        event%icu, &
        event%izone, &
        event%istatus, &
        event%get_code(), &
        event%trelease, &
        event%ttrack, &
        event%x, &
        event%y, &
        event%z, &
        particle%name
    end if
  end subroutine save_event

  !> @brief Is the output unit valid?
  logical function should_print(this)
    class(ParticleTracksType), intent(inout) :: this
    should_print = this%iout >= 0
  end function should_print

  !> @brief Write a particle event to files for which the particle
  !! is eligible, and print the event to output unit if requested.
  !! This function is the module's main entry point. It should be
  !! subscribed as an event handler to particle event dispatchers.
  function write_particle_event(context, particle, event) result(handled)
    ! dummy
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: handled
    ! local
    integer(I4B) :: i
    type(ParticleTrackFileType) :: file

    select type (context)
    type is (ParticleTracksType)
      if (context%should_print()) &
        call event%log(context%iout)
      if (context%is_selected(event)) then
        do i = 1, context%ntrackfiles
          file = context%files(i)
          if (context%should_save(particle, file)) &
            call save_event(file%iun, particle, event, csv=file%csv)
        end do
      end if
      handled = .true.
    end select
  end function write_particle_event

end module ParticleTracksModule

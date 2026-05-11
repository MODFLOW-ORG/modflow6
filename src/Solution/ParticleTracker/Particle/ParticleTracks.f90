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
!!
!! Events can be buffered in memory or in a scratch file, and in either
!! case are flushed to disk only when a time step successfully finishes.
!<
module ParticleTracksModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType, ACTIVE
  use ParticleEventModule, only: ParticleEventType
  use ReleaseEventModule, only: ReleaseEventType
  use TimeStepEventModule, only: TimeStepEventType
  use TerminationEventModule, only: TerminationEventType
  use WeakSinkEventModule, only: WeakSinkEventType
  use UserTimeEventModule, only: UserTimeEventType
  use CellExitEventModule, only: CellExitEventType
  use SubcellExitEventModule, only: SubcellExitEventType
  use DroppedEventModule, only: DroppedEventType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform
  use ParticleTrackEventBufferModule
  use MemoryBufferModule, only: MemoryBufferType
  use ScratchFileBufferModule, only: ScratchFileBufferType

  implicit none
  public :: ParticleTracksType, &
            ParticleTrackEventSelectionType, &
            add_particle_event

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
  !! that only certain event types are printed or written to files. Particle
  !! events are buffered in memory or in a scratch file, flushed to disk only
  !! when the time step is successfully solved for the last time (there may
  !! be multiple solves per time step, depending on ATS and Picard options).
  !<
  type :: ParticleTracksType
    private
    integer(I4B), public :: iout = -1 !< log file unit
    integer(I4B), public :: ntrackfiles !< number of track files
    type(ParticleTrackFileType), public, allocatable :: files(:) !< track files
    type(ParticleTrackEventSelectionType), public :: selected !< event selection
    class(ParticleTrackEventBufferType), allocatable :: buffer !< event buffer
  contains
    procedure, public :: init_file
    procedure, public :: init_buffer
    procedure, public :: is_selected
    procedure, public :: select_events
    procedure, public :: buffer_event
    procedure, public :: flush_buffer
    procedure, public :: discard_buffer
    procedure, public :: destroy
    procedure :: expand_files
  end type ParticleTracksType

contains

  !> @brief Initialize a binary or CSV output file.
  subroutine init_file(this, iun, csv, iprp)
    class(ParticleTracksType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    type(ParticleTrackFileType), pointer :: file

    call this%expand_files(increment=1)

    allocate (file)
    file%iun = iun
    if (present(csv)) file%csv = csv
    if (present(iprp)) file%iprp = iprp
    this%ntrackfiles = size(this%files)
    this%files(this%ntrackfiles) = file
  end subroutine init_file

  !> @brief Initialize the event buffer strategy
  subroutine init_buffer(this, scratch)
    class(ParticleTracksType) :: this
    logical(LGP), intent(in) :: scratch

    if (scratch) then
      allocate (ScratchFileBufferType :: this%buffer)
    else
      allocate (MemoryBufferType :: this%buffer)
    end if
    call this%buffer%init()
  end subroutine init_buffer

  !> @brief Destroy the particle track manager.
  subroutine destroy(this)
    class(ParticleTracksType) :: this
    call this%buffer%destroy()
    if (allocated(this%buffer)) deallocate (this%buffer)
    if (allocated(this%files)) deallocate (this%files)
  end subroutine destroy

  !> @brief Grow the array of track files.
  subroutine expand_files(this, increment)
    class(ParticleTracksType) :: this
    integer(I4B), optional, intent(in) :: increment
    integer(I4B) :: inclocal, isize, newsize
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

  !> @brief Buffer an event for deferred write.
  subroutine buffer_event(this, particle, event)
    class(ParticleTracksType) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    type(TrackRecordType) :: rec

    rec%kper = event%kper; rec%kstp = event%kstp
    rec%imdl = event%imdl; rec%iprp = event%iprp
    rec%irpt = event%irpt; rec%ilay = event%ilay
    rec%icu = event%icu; rec%izone = event%izone
    rec%istatus = event%istatus; rec%ireason = event%get_code()
    rec%trelease = event%trelease; rec%ttrack = event%ttrack
    rec%x = event%x; rec%y = event%y; rec%z = event%z
    rec%name = trim(adjustl(particle%name))

    call this%buffer%append(rec)
  end subroutine buffer_event

  !> @brief Flush the event buffer to disk.
  subroutine flush_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%flush(this%files)
  end subroutine flush_buffer

  !> @brief Discard buffered events without writing.
  subroutine discard_buffer(this)
    class(ParticleTracksType) :: this
    call this%buffer%discard()
  end subroutine discard_buffer

  !> @brief Add a particle event to be written to eligible
  !! files and printed to an output file unit if requested.
  !! This function should be subscribed as an event handler
  !! to particle event dispatchers. Events are buffered, to
  !! be written to output files upon successful completion
  !! of a time step when the framework OT hook is executed.
  function add_particle_event(context, particle, event) result(handled)
    class(*), pointer :: context
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer, intent(in) :: event
    logical(LGP) :: handled

    select type (context)
    type is (ParticleTracksType)
      if (context%iout >= 0) &
        call event%log(context%iout)
      if (context%is_selected(event)) &
        call context%buffer_event(particle, event)
      handled = .true.
    end select
  end function add_particle_event

end module ParticleTracksModule

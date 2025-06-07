module TrackControlModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, DPIO180
  use ParticleModule, only: ParticleType
  use BaseDisModule, only: DisBaseType
  use GeomUtilModule, only: transform
  use TrackFileModule, only: TrackFileType, save_record

  implicit none
  public :: TrackControlType

  !> @brief Manages particle track (i.e. pathline) files.
  !!
  !! Optionally filters events ("ireason" codes, selectable in the PRT-OC pkg):
  !!
  !!    0: RELEASE: particle is released
  !!    1: TRANSIT: particle moves from cell to cell
  !!    2: TIMESTEP: timestep ends
  !!    3: TERMINATE: tracking stops for a particle
  !!    4: WEAKSINK: particle exits a weak sink
  !!    5: USERTIME: user-specified tracking time
  !!
  !! An arbitrary number of files can be managed. Internal arrays
  !! are resized as needed.
  !<
  type :: TrackControlType
    private
    type(TrackFileType), public, allocatable :: trackfiles(:) !< output files
    integer(I4B), public :: ntrackfiles !< number of output files
    logical(LGP), public :: trackrelease !< track release events
    logical(LGP), public :: trackexit !< track cell-to-cell transitions
    logical(LGP), public :: tracktimestep !< track timestep ends
    logical(LGP), public :: trackterminate !< track termination events
    logical(LGP), public :: trackweaksink !< track weak sink exit events
    logical(LGP), public :: trackusertime !< track user-selected times
  contains
    procedure :: expand
    procedure, public :: init_track_file
    procedure, public :: save
    procedure, public :: set_track_events
  end type TrackControlType

contains

  !> @brief Initialize a new track file
  subroutine init_track_file(this, iun, csv, iprp)
    ! dummy
    class(TrackControlType) :: this
    integer(I4B), intent(in) :: iun
    logical(LGP), intent(in), optional :: csv
    integer(I4B), intent(in), optional :: iprp
    ! local
    type(TrackFileType), pointer :: file

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

  !> @brief Expand the trackfile array, internal use only
  subroutine expand(this, increment)
    ! dummy
    class(TrackControlType) :: this
    integer(I4B), optional, intent(in) :: increment
    ! local
    integer(I4B) :: inclocal
    integer(I4B) :: isize
    integer(I4B) :: newsize
    type(TrackFileType), allocatable, dimension(:) :: temp

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

  !> @brief Save the particle's state to track output file(s).
  !!
  !! A record is saved to all enabled model-level files and to
  !! any PRP-level files with PRP index matching the particle's
  !! PRP index.
  !<
  subroutine save(this, particle, dis, kper, kstp, reason, level)
    ! dummy
    class(TrackControlType), intent(inout) :: this
    type(ParticleType), pointer, intent(in) :: particle
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: kper
    integer(I4B), intent(in) :: kstp
    integer(I4B), intent(in) :: reason
    integer(I4B), intent(in), optional :: level
    ! local
    integer(I4B) :: i
    type(TrackFileType) :: file

    ! Only save if reporting is enabled for specified event.
    if (.not. ((this%trackrelease .and. reason == 0) .or. &
               (this%trackexit .and. reason == 1) .or. &
               (this%tracktimestep .and. reason == 2) .or. &
               (this%trackterminate .and. reason == 3) .or. &
               (this%trackweaksink .and. reason == 4) .or. &
               (this%trackusertime .and. reason == 5))) &
      return

    ! For now, only allow reporting from outside the tracking
    ! algorithm (e.g. release time), in which case level will
    ! not be provided, or if within the tracking solution, in
    ! subcells (level 3) only. This may change if the subcell
    ! ever delegates tracking to even smaller subcomponents.
    if (present(level)) then
      if (level .ne. 3) return
    end if

    ! Save to any enabled model-scoped or PRP-scoped files
    do i = 1, this%ntrackfiles
      file = this%trackfiles(i)
      if (file%iun > 0 .and. &
          (file%iprp == -1 .or. &
           file%iprp == particle%iprp)) &
        call save_record(file%iun, particle, dis, &
                         kper, kstp, reason, csv=file%csv)
    end do
  end subroutine save

  !> @brief Configure particle events to track.
  !!
  !! Each tracking event corresponds to an "ireason" code
  !! as appears in each row of track output.
  !<
  subroutine set_track_events(this, &
                              release, &
                              cellexit, &
                              timestep, &
                              terminate, &
                              weaksink, &
                              usertime)
    class(TrackControlType) :: this
    logical(LGP), intent(in) :: release
    logical(LGP), intent(in) :: cellexit
    logical(LGP), intent(in) :: timestep
    logical(LGP), intent(in) :: terminate
    logical(LGP), intent(in) :: weaksink
    logical(LGP), intent(in) :: usertime
    this%trackrelease = release
    this%trackexit = cellexit
    this%tracktimestep = timestep
    this%trackterminate = terminate
    this%trackweaksink = weaksink
    this%trackusertime = usertime
  end subroutine set_track_events

end module TrackControlModule

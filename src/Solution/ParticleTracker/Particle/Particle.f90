module ParticleModule

  use KindModule, only: DP, I4B, LGP
  use ListModule, only: ListType
  use ConstantsModule, only: DZERO, DONE, LENMEMPATH, LENBOUNDNAME, &
                             LINELENGTH
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_reallocate
  implicit none
  public

  !> Tracking "levels" defined in method modules. Currently only 3 used.
  integer, parameter :: MAX_LEVEL = 4

  !> @brief Particle status enumeration.
  !!
  !! Particles begin in status 1 (active) at release time. Status may only
  !! increase over time.  Status values greater than one imply termination.
  !! A particle may terminate for several reasons, all mutually exclusive.
  !! A particle's final tracking status will always be greater than one.
  !!
  !! Status codes 0-3 and 5-8 correspond directly to MODPATH 7 status codes.
  !! Code 4 does not apply to PRT because PRT does not distinguish forwards
  !! from backwards tracking. Status code 9 provides more specific, subcell-
  !! level information about a particle which terminates due to no outflow.
  !! Code 10 distinguishes particles which have "timed out" upon reaching a
  !! user-specified stop time or the end of the simulation.
  !<
  enum, bind(C)
    enumerator :: ACTIVE = 1
    enumerator :: TERM_BOUNDARY = 2 !< terminated at a boundary face
    enumerator :: TERM_WEAKSINK = 3 !< terminated in a weak sink cell
    enumerator :: TERM_NO_EXITS = 5 !< terminated in a cell with no exit face
    enumerator :: TERM_STOPZONE = 6 !< terminated in a cell with a stop zone number
    enumerator :: TERM_INACTIVE = 7 !< terminated in an inactive cell
    enumerator :: TERM_UNRELEASED = 8 !< terminated permanently unreleased
    enumerator :: TERM_NO_EXITS_SUB = 9 !< terminated in a subcell with no exit face
    enumerator :: TERM_TIMEOUT = 10 !< terminated at stop time or end of simulation
  end enum

  !> @brief Particle tracked by the PRT model.
  !!
  !! Record-type to conveniently shuffle a particle's
  !! state to/from storage before/after its trajectory
  !! is solved for each time step.
  !!
  !! Particle coordinates may be local to the cell or
  !! global/model. Routines are provided to convert a
  !! particle's global coordinates to/from cell-local
  !! coordinates for tracking through cell subdomains.
  !!
  !! Particles are identified by composite key, i.e.,
  !! combinations of properties imdl, iprp, irpt, and
  !! trelease. An optional label may be provided, but
  !! need not be unique.
  !!
  !! Particles carry a coordinate transform along with
  !! themselves, instead of just letting the tracking
  !! methods handle transforms, because we may need to
  !! report an event at any time in model coordinates.
  !<
  type ParticleType
    private
    ! identity
    character(len=LENBOUNDNAME), public :: name = '' !< optional particle name
    integer(I4B), public :: imdl !< index of model the particle originated in
    integer(I4B), public :: iprp !< index of release package the particle is from
    integer(I4B), public :: irpt !< index of release point the particle is from
    integer(I4B), public :: ip !< index of particle in the particle list
    ! options
    logical(LGP), public :: extend !< whether to extend tracking beyond the end of the simulation
    logical(LGP), public :: frctrn !< whether to force solving the particle with the ternary method
    integer(I4B), public :: istopweaksink !< weak sink option (0: do not stop, 1: stop)
    integer(I4B), public :: istopzone !< stop zone number
    integer(I4B), public :: idrymeth !< dry tracking method
    integer(I4B), public :: iexmeth !< method for iterative solution of particle exit location and time in generalized Pollock's method
    integer(I4B), public :: icycwin !< cycle detection window size
    real(DP), public :: extol !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    ! state
    integer(I4B), public :: itrdomain(MAX_LEVEL) !< tracking domain indices
    integer(I4B), public :: iboundary(MAX_LEVEL) !< tracking domain boundary indices
    integer(I4B), public :: icu !< user cell number
    integer(I4B), public :: ilay !< grid layer
    integer(I4B), public :: izone !< current zone number
    integer(I4B), public :: istatus !< tracking status
    real(DP), public :: x !< x coordinate
    real(DP), public :: y !< y coordinate
    real(DP), public :: z !< z coordinate
    real(DP), public :: trelease !< release time
    real(DP), public :: tstop !< stop time
    real(DP), public :: ttrack !< time tracked so far
    real(DP), public :: xorigin !< x origin for coordinate transformation from model to local
    real(DP), public :: yorigin !< y origin for coordinate transformation from model to local
    real(DP), public :: zorigin !< z origin for coordinate transformation from model to local
    real(DP), public :: sinrot !< sine of rotation angle for coordinate transformation from model to local
    real(DP), public :: cosrot !< cosine of rotation angle for coordinate transformation from model to local
    logical(LGP), public :: transformed !< whether coordinates have been transformed from model to local
    logical(LGP), public :: advancing !< whether particle is still being tracked for current time step
    type(ListType), public, pointer :: history !< history of particle positions (for cycle detection)
  contains
    procedure, public :: destroy => destroy_particle
    procedure, public :: get_model_coords
    procedure, public :: transform => transform_coords
    procedure, public :: reset_transform
    procedure, public :: get_id
  end type ParticleType

  !> @brief Structure of arrays to store particles.
  type ParticleStoreType
    private
    ! identity
    character(len=LENBOUNDNAME), dimension(:), pointer, public, contiguous :: name !< optional particle label
    integer(I4B), dimension(:), pointer, public, contiguous :: imdl !< index of model particle originated in
    integer(I4B), dimension(:), pointer, public, contiguous :: iprp !< index of release package the particle originated in
    integer(I4B), dimension(:), pointer, public, contiguous :: irpt !< index of release point in the particle release package the particle originated in
    ! options
    logical(LGP), dimension(:), pointer, public, contiguous :: extend !< whether to extend tracking beyond the end of the simulation
    logical(LGP), dimension(:), pointer, public, contiguous :: frctrn !< force ternary method
    integer(I4B), dimension(:), pointer, public, contiguous :: istopweaksink !< weak sink option: 0 = do not stop, 1 = stop
    integer(I4B), dimension(:), pointer, public, contiguous :: istopzone !< stop zone number
    integer(I4B), dimension(:), pointer, public, contiguous :: idrymeth !< stop in dry cells
    integer(I4B), dimension(:), pointer, public, contiguous :: iexmeth !< method for iterative solution of particle exit location and time in generalized Pollock's method
    integer(I4B), dimension(:), pointer, public, contiguous :: icycwin !< cycle detection window size
    real(DP), dimension(:), pointer, public, contiguous :: extol !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    ! state
    integer(I4B), dimension(:, :), pointer, public, contiguous :: itrdomain !< array of indices for domains in the tracking domain hierarchy
    integer(I4B), dimension(:, :), pointer, public, contiguous :: iboundary !< array of indices for tracking domain boundaries
    integer(I4B), dimension(:), pointer, public, contiguous :: icu !< cell number (user)
    integer(I4B), dimension(:), pointer, public, contiguous :: ilay !< layer
    integer(I4B), dimension(:), pointer, public, contiguous :: izone !< current zone number
    integer(I4B), dimension(:), pointer, public, contiguous :: istatus !< particle status
    real(DP), dimension(:), pointer, public, contiguous :: x !< model x coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: y !< model y coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: z !< model z coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: trelease !< particle release time
    real(DP), dimension(:), pointer, public, contiguous :: tstop !< particle stop time
    real(DP), dimension(:), pointer, public, contiguous :: ttrack !< current tracking time
  contains
    procedure, public :: destroy
    procedure, public :: num_stored
    procedure, public :: resize
    procedure, public :: get
    procedure, public :: put
  end type ParticleStoreType

contains

  !> @brief Create a new particle
  subroutine create_particle(particle)
    type(ParticleType), pointer :: particle !< particle
    allocate (particle)
    allocate (particle%history)
  end subroutine create_particle

  !> @brief Allocate particle store
  subroutine create_particle_store(store, np, mempath)
    type(ParticleStoreType), pointer :: store !< store
    integer(I4B), intent(in) :: np !< number of particles
    character(*), intent(in) :: mempath !< path to memory

    allocate (store)
    call mem_allocate(store%imdl, np, 'PLIMDL', mempath)
    call mem_allocate(store%irpt, np, 'PLIRPT', mempath)
    call mem_allocate(store%iprp, np, 'PLIPRP', mempath)
    call mem_allocate(store%name, LENBOUNDNAME, np, 'PLNAME', mempath)
    call mem_allocate(store%icu, np, 'PLICU', mempath)
    call mem_allocate(store%ilay, np, 'PLILAY', mempath)
    call mem_allocate(store%izone, np, 'PLIZONE', mempath)
    call mem_allocate(store%istatus, np, 'PLISTATUS', mempath)
    call mem_allocate(store%x, np, 'PLX', mempath)
    call mem_allocate(store%y, np, 'PLY', mempath)
    call mem_allocate(store%z, np, 'PLZ', mempath)
    call mem_allocate(store%trelease, np, 'PLTRELEASE', mempath)
    call mem_allocate(store%tstop, np, 'PLTSTOP', mempath)
    call mem_allocate(store%ttrack, np, 'PLTTRACK', mempath)
    call mem_allocate(store%istopweaksink, np, 'PLISTOPWEAKSINK', mempath)
    call mem_allocate(store%istopzone, np, 'PLISTOPZONE', mempath)
    call mem_allocate(store%idrymeth, np, 'PLIDRYMETH', mempath)
    call mem_allocate(store%frctrn, np, 'PLFRCTRN', mempath)
    call mem_allocate(store%iexmeth, np, 'PLIEXMETH', mempath)
    call mem_allocate(store%extol, np, 'PLEXTOL', mempath)
    call mem_allocate(store%extend, np, 'PLEXTEND', mempath)
    call mem_allocate(store%icycwin, np, 'PLICYCWIN', mempath)
    call mem_allocate(store%itrdomain, np, MAX_LEVEL, 'PLIDOMAIN', mempath)
    call mem_allocate(store%iboundary, np, MAX_LEVEL, 'PLIBOUNDARY', mempath)
  end subroutine create_particle_store

  !> @brief Destroy particle store after use.
  subroutine destroy(this, mempath)
    class(ParticleStoreType), intent(inout) :: this !< store
    character(*), intent(in) :: mempath !< path to memory

    call mem_deallocate(this%imdl, 'PLIMDL', mempath)
    call mem_deallocate(this%iprp, 'PLIPRP', mempath)
    call mem_deallocate(this%irpt, 'PLIRPT', mempath)
    call mem_deallocate(this%name, 'PLNAME', mempath)
    call mem_deallocate(this%icu, 'PLICU', mempath)
    call mem_deallocate(this%ilay, 'PLILAY', mempath)
    call mem_deallocate(this%izone, 'PLIZONE', mempath)
    call mem_deallocate(this%istatus, 'PLISTATUS', mempath)
    call mem_deallocate(this%x, 'PLX', mempath)
    call mem_deallocate(this%y, 'PLY', mempath)
    call mem_deallocate(this%z, 'PLZ', mempath)
    call mem_deallocate(this%trelease, 'PLTRELEASE', mempath)
    call mem_deallocate(this%tstop, 'PLTSTOP', mempath)
    call mem_deallocate(this%ttrack, 'PLTTRACK', mempath)
    call mem_deallocate(this%istopweaksink, 'PLISTOPWEAKSINK', mempath)
    call mem_deallocate(this%istopzone, 'PLISTOPZONE', mempath)
    call mem_deallocate(this%idrymeth, 'PLIDRYMETH', mempath)
    call mem_deallocate(this%frctrn, 'PLFRCTRN', mempath)
    call mem_deallocate(this%iexmeth, 'PLIEXMETH', mempath)
    call mem_deallocate(this%extol, 'PLEXTOL', mempath)
    call mem_deallocate(this%extend, 'PLEXTEND', mempath)
    call mem_deallocate(this%icycwin, 'PLICYCWIN', mempath)
    call mem_deallocate(this%itrdomain, 'PLIDOMAIN', mempath)
    call mem_deallocate(this%iboundary, 'PLIBOUNDARY', mempath)
  end subroutine destroy

  !> @brief Destroy a particle after use.
  subroutine destroy_particle(particle)
    class(ParticleType), intent(inout) :: particle !< particle
    deallocate (particle%history)
  end subroutine destroy_particle

  !> @brief Reallocate particle storage to the given size.
  subroutine resize(this, np, mempath)
    ! dummy
    class(ParticleStoreType), intent(inout) :: this !< particle store
    integer(I4B), intent(in) :: np !< number of particles
    character(*), intent(in) :: mempath !< path to memory

    ! resize arrays
    call mem_reallocate(this%imdl, np, 'PLIMDL', mempath)
    call mem_reallocate(this%iprp, np, 'PLIPRP', mempath)
    call mem_reallocate(this%irpt, np, 'PLIRPT', mempath)
    call mem_reallocate(this%name, LENBOUNDNAME, np, 'PLNAME', mempath)
    call mem_reallocate(this%icu, np, 'PLICU', mempath)
    call mem_reallocate(this%ilay, np, 'PLILAY', mempath)
    call mem_reallocate(this%izone, np, 'PLIZONE', mempath)
    call mem_reallocate(this%istatus, np, 'PLISTATUS', mempath)
    call mem_reallocate(this%x, np, 'PLX', mempath)
    call mem_reallocate(this%y, np, 'PLY', mempath)
    call mem_reallocate(this%z, np, 'PLZ', mempath)
    call mem_reallocate(this%trelease, np, 'PLTRELEASE', mempath)
    call mem_reallocate(this%tstop, np, 'PLTSTOP', mempath)
    call mem_reallocate(this%ttrack, np, 'PLTTRACK', mempath)
    call mem_reallocate(this%istopweaksink, np, 'PLISTOPWEAKSINK', mempath)
    call mem_reallocate(this%istopzone, np, 'PLISTOPZONE', mempath)
    call mem_reallocate(this%idrymeth, np, 'PLIDRYMETH', mempath)
    call mem_reallocate(this%frctrn, np, 'PLFRCTRN', mempath)
    call mem_reallocate(this%iexmeth, np, 'PLIEXMETH', mempath)
    call mem_reallocate(this%extol, np, 'PLEXTOL', mempath)
    call mem_reallocate(this%extend, np, 'PLEXTEND', mempath)
    call mem_reallocate(this%icycwin, np, 'PLICYCWIN', mempath)
    call mem_reallocate(this%itrdomain, np, MAX_LEVEL, 'PLIDOMAIN', mempath)
    call mem_reallocate(this%iboundary, np, MAX_LEVEL, 'PLIBOUNDARY', mempath)
  end subroutine resize

  !> @brief Load a particle from the particle store.
  !!
  !! This routine is used to initialize a particle for tracking.
  !! The advancing flag and coordinate transformation are reset.
  !<
  subroutine get(this, particle, imdl, iprp, ip)
    class(ParticleStoreType), intent(inout) :: this !< particle store
    class(ParticleType), intent(inout) :: particle !< particle
    integer(I4B), intent(in) :: imdl !< index of model particle originated in
    integer(I4B), intent(in) :: iprp !< index of particle release package particle originated in
    integer(I4B), intent(in) :: ip !< index into the particle list

    call particle%reset_transform()
    call particle%history%Clear()
    particle%imdl = imdl
    particle%iprp = iprp
    particle%irpt = this%irpt(ip)
    particle%ip = ip
    particle%name = this%name(ip)
    particle%istopweaksink = this%istopweaksink(ip)
    particle%istopzone = this%istopzone(ip)
    particle%idrymeth = this%idrymeth(ip)
    particle%icu = this%icu(ip)
    particle%ilay = this%ilay(ip)
    particle%izone = this%izone(ip)
    particle%istatus = this%istatus(ip)
    particle%x = this%x(ip)
    particle%y = this%y(ip)
    particle%z = this%z(ip)
    particle%trelease = this%trelease(ip)
    particle%tstop = this%tstop(ip)
    particle%ttrack = this%ttrack(ip)
    particle%advancing = .true.
    particle%itrdomain(1:MAX_LEVEL) = &
      this%itrdomain(ip, 1:MAX_LEVEL)
    particle%itrdomain(1) = imdl
    particle%iboundary(1:MAX_LEVEL) = &
      this%iboundary(ip, 1:MAX_LEVEL)
    particle%frctrn = this%frctrn(ip)
    particle%iexmeth = this%iexmeth(ip)
    particle%extol = this%extol(ip)
    particle%extend = this%extend(ip)
    particle%icycwin = this%icycwin(ip)
  end subroutine get

  !> @brief Save a particle's state to the particle store.
  subroutine put(this, particle, ip)
    class(ParticleStoreType), intent(inout) :: this !< particle storage
    class(ParticleType), intent(in) :: particle !< particle
    integer(I4B), intent(in) :: ip !< particle index

    this%imdl(ip) = particle%imdl
    this%iprp(ip) = particle%iprp
    this%irpt(ip) = particle%irpt
    this%name(ip) = particle%name
    this%istopweaksink(ip) = particle%istopweaksink
    this%istopzone(ip) = particle%istopzone
    this%idrymeth(ip) = particle%idrymeth
    this%icu(ip) = particle%icu
    this%ilay(ip) = particle%ilay
    this%izone(ip) = particle%izone
    this%istatus(ip) = particle%istatus
    this%x(ip) = particle%x
    this%y(ip) = particle%y
    this%z(ip) = particle%z
    this%trelease(ip) = particle%trelease
    this%tstop(ip) = particle%tstop
    this%ttrack(ip) = particle%ttrack
    this%itrdomain( &
      ip, &
      1:MAX_LEVEL) = &
      particle%itrdomain(1:MAX_LEVEL)
    this%iboundary( &
      ip, &
      1:MAX_LEVEL) = &
      particle%iboundary(1:MAX_LEVEL)
    this%frctrn(ip) = particle%frctrn
    this%iexmeth(ip) = particle%iexmeth
    this%extol(ip) = particle%extol
    this%extend(ip) = particle%extend
    this%icycwin(ip) = particle%icycwin
  end subroutine put

  !> @brief Transform particle coordinates.
  !!
  !! Apply a translation and/or rotation to particle coordinates.
  !! No rescaling. It's also possible to invert a transformation.
  !! Be sure to reset the transformation after using it.
  !<
  subroutine transform_coords(this, xorigin, yorigin, zorigin, &
                              sinrot, cosrot, invert)
    use GeomUtilModule, only: transform, compose
    class(ParticleType), intent(inout) :: this !< particle
    real(DP), intent(in), optional :: xorigin !< x coordinate of origin
    real(DP), intent(in), optional :: yorigin !< y coordinate of origin
    real(DP), intent(in), optional :: zorigin !< z coordinate of origin
    real(DP), intent(in), optional :: sinrot !< sine of rotation angle
    real(DP), intent(in), optional :: cosrot !< cosine of rotation angle
    logical(LGP), intent(in), optional :: invert !< whether to invert

    call transform(this%x, this%y, this%z, &
                   this%x, this%y, this%z, &
                   xorigin, yorigin, zorigin, &
                   sinrot, cosrot, invert)

    ! Compose is needed only because we have to untransform
    ! coordinates: we may need to report a particle event
    ! at any point in the tracking method hierarchy, so we
    ! need to know how to map the particle position back to
    ! model coords from coords local to the current domain.
    call compose(this%xorigin, this%yorigin, this%zorigin, &
                 this%sinrot, this%cosrot, &
                 xorigin, yorigin, zorigin, &
                 sinrot, cosrot, invert)

    this%transformed = .true.
  end subroutine transform_coords

  !> @brief Reset particle coordinate transformation properties.
  subroutine reset_transform(this)
    class(ParticleType), intent(inout) :: this !< particle

    this%xorigin = DZERO
    this%yorigin = DZERO
    this%zorigin = DZERO
    this%sinrot = DZERO
    this%cosrot = DONE
    this%cosrot = DONE
    this%transformed = .false.
  end subroutine reset_transform

  !> @brief Return the particle's model coordinates,
  !! inverting any applied transformation if needed.
  !! The particle's state is not altered.
  subroutine get_model_coords(this, x, y, z)
    use GeomUtilModule, only: transform, compose
    class(ParticleType), intent(in) :: this !< particle
    real(DP), intent(out) :: x !< x coordinate
    real(DP), intent(out) :: y !< y coordinate
    real(DP), intent(out) :: z !< z coordinate

    if (this%transformed) then
      call transform(this%x, this%y, this%z, x, y, z, &
                     this%xorigin, this%yorigin, this%zorigin, &
                     this%sinrot, this%cosrot, invert=.true.)
    else
      x = this%x
      y = this%y
      z = this%z
    end if
  end subroutine get_model_coords

  !> @brief Return the number of particles.
  integer function num_stored(this) result(n)
    class(ParticleStoreType) :: this
    n = size(this%imdl)
  end function num_stored

  !> @brief Get a string identifier for the particle.
  function get_id(this) result(str)
    class(ParticleType), intent(in) :: this
    character(len=:), allocatable :: str
    ! local
    character(len=LINELENGTH) :: temp

    write (temp, '(I0,1a,I0,1a,I0,1a,G0)') &
      this%imdl, this%iprp, this%irpt, this%trelease
    str = trim(adjustl(temp))
  end function get_id

end module ParticleModule

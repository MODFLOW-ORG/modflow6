module ParticleModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LENMEMPATH, LENBOUNDNAME
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_reallocate
  implicit none
  public

  !> Tracking "levels" (1: model, 2: cell, 3: subcell). A
  !! level identifies the domain through which a tracking
  !! method is responsible for moving a particle. Methods
  !! each operate on a particular level, delegating among
  !! more methods as appropriate for finer-grained levels.
  integer, parameter :: MAX_LEVEL = 4

  !> @brief Particle tracking level.
  !!
  !! Identifies the domain over which a tracking method is responsible for
  !! moving a particle. Tracking methods at lower levels (coarser-grained)
  !! defer to higher-level (finer-grained) methods within their subdomains.
  !<
  enum, bind(C)
    enumerator :: LVL_MODEL = 1
    enumerator :: LVL_CELL = 2
    enumerator :: LVL_SUBCELL = 3
  end enum

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

  !> @brief Particle event enumeration.
  !!
  !! A number of events may occur to particles, each of which may (or may
  !! not) be of interest to the user. The user selects among events to be
  !! reported. A corresponding event code is reported with each record to
  !! identify the record's cause.
  !!
  !! Records may be identical except for their event code, reflecting the
  !! fact that multiple events of interest may occur at any given moment.
  !<
  enum, bind(C)
    enumerator :: RELEASE = 0 !< particle was released
    enumerator :: EXIT = 1 !< particle exited a cell
    enumerator :: TIMESTEP = 2 !< time step ended
    enumerator :: TERMINATE = 3 !< particle terminated
    enumerator :: WEAKSINK = 4 !< particle entered a weak sink cell
    enumerator :: USERTIME = 5 !< user-specified tracking time
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
  !! need not be unique
  !<
  type ParticleType
    private
    ! identity
    character(len=LENBOUNDNAME), public :: name = '' !< optional particle name
    integer(I4B), public :: imdl !< index of model the particle originated in
    integer(I4B), public :: iprp !< index of release package the particle is from
    integer(I4B), public :: irpt !< index of release point the particle is from
    integer(I4B), public :: ip !< index of particle in the particle list
    ! stop criteria
    integer(I4B), public :: istopweaksink !< weak sink option (0: do not stop, 1: stop)
    integer(I4B), public :: istopzone !< stop zone number
    integer(I4B), public :: idrymeth !< dry tracking method
    ! state
    integer(I4B), allocatable, public :: itrdomain(:) !< tracking domain
    integer(I4B), allocatable, public :: iboundary(:) !< tracking domain boundaries
    integer(I4B), public :: icp !< previous cell number (reduced)
    integer(I4B), public :: icu !< user cell number
    integer(I4B), public :: ilay !< grid layer
    integer(I4B), public :: izone !< current zone number
    integer(I4B), public :: izp !< previous zone number
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
    real(DP), public :: extol !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    logical(LGP), public :: transformed !< whether coordinates have been transformed from model to local
    logical(LGP), public :: advancing !< whether particle is still being tracked for current time step
    integer(I4B), public :: ifrctrn !< whether to force solving the particle with the ternary method
    integer(I4B), public :: iexmeth !< method for iterative solution of particle exit location and time in generalized Pollock's method
    integer(I4B), public :: iextend !< whether to extend tracking beyond the end of the simulation
  contains
    procedure, public :: get_model_coords
    procedure, public :: transform => transform_coords
    procedure, public :: reset_transform
  end type ParticleType

  !> @brief Structure of arrays to store particles.
  type ParticleStoreType
    private
    ! identity
    character(len=LENBOUNDNAME), dimension(:), pointer, public, contiguous :: name !< optional particle label
    integer(I4B), dimension(:), pointer, public, contiguous :: imdl !< index of model particle originated in
    integer(I4B), dimension(:), pointer, public, contiguous :: iprp !< index of release package the particle originated in
    integer(I4B), dimension(:), pointer, public, contiguous :: irpt !< index of release point in the particle release package the particle originated in
    ! stopping criteria
    integer(I4B), dimension(:), pointer, public, contiguous :: istopweaksink !< weak sink option: 0 = do not stop, 1 = stop
    integer(I4B), dimension(:), pointer, public, contiguous :: istopzone !< stop zone number
    integer(I4B), dimension(:), pointer, public, contiguous :: idrymeth !< stop in dry cells
    ! state
    integer(I4B), dimension(:, :), pointer, public, contiguous :: itrdomain !< array of indices for domains in the tracking domain hierarchy
    integer(I4B), dimension(:, :), pointer, public, contiguous :: iboundary !< array of indices for tracking domain boundaries
    integer(I4B), dimension(:), pointer, public, contiguous :: icu !< cell number (user)
    integer(I4B), dimension(:), pointer, public, contiguous :: ilay !< layer
    integer(I4B), dimension(:), pointer, public, contiguous :: izone !< current zone number
    integer(I4B), dimension(:), pointer, public, contiguous :: izp !< previous zone number
    integer(I4B), dimension(:), pointer, public, contiguous :: istatus !< particle status
    real(DP), dimension(:), pointer, public, contiguous :: x !< model x coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: y !< model y coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: z !< model z coord of particle
    real(DP), dimension(:), pointer, public, contiguous :: trelease !< particle release time
    real(DP), dimension(:), pointer, public, contiguous :: tstop !< particle stop time
    real(DP), dimension(:), pointer, public, contiguous :: ttrack !< current tracking time
    integer(I4B), dimension(:), pointer, public, contiguous :: ifrctrn !< force ternary method
    integer(I4B), dimension(:), pointer, public, contiguous :: iexmeth !< method for iterative solution of particle exit location and time in generalized Pollock's method
    real(DP), dimension(:), pointer, public, contiguous :: extol !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    integer(LGP), dimension(:), pointer, public, contiguous :: extend !< whether to extend tracking beyond the end of the simulation
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
    allocate (particle%itrdomain(MAX_LEVEL))
    allocate (particle%iboundary(MAX_LEVEL))
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
    call mem_allocate(store%izp, np, 'PLIZP', mempath)
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
    call mem_allocate(store%ifrctrn, np, 'PLIFRCTRN', mempath)
    call mem_allocate(store%iexmeth, np, 'PLIEXMETH', mempath)
    call mem_allocate(store%extol, np, 'PLEXTOL', mempath)
    call mem_allocate(store%extend, np, 'PLIEXTEND', mempath)
    call mem_allocate(store%itrdomain, np, MAX_LEVEL, 'PLITRDOMAIN', mempath)
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
    call mem_deallocate(this%izp, 'PLIZP', mempath)
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
    call mem_deallocate(this%ifrctrn, 'PLIFRCTRN', mempath)
    call mem_deallocate(this%iexmeth, 'PLIEXMETH', mempath)
    call mem_deallocate(this%extol, 'PLEXTOL', mempath)
    call mem_deallocate(this%extend, 'PLIEXTEND', mempath)
    call mem_deallocate(this%itrdomain, 'PLITRDOMAIN', mempath)
    call mem_deallocate(this%iboundary, 'PLIBOUNDARY', mempath)
  end subroutine destroy

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
    call mem_reallocate(this%izp, np, 'PLIZP', mempath)
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
    call mem_reallocate(this%ifrctrn, np, 'PLIFRCTRN', mempath)
    call mem_reallocate(this%iexmeth, np, 'PLIEXMETH', mempath)
    call mem_reallocate(this%extol, np, 'PLEXTOL', mempath)
    call mem_reallocate(this%extend, np, 'PLIEXTEND', mempath)
    call mem_reallocate(this%itrdomain, np, MAX_LEVEL, 'PLITRDOMAIN', mempath)
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
    particle%imdl = imdl
    particle%iprp = iprp
    particle%irpt = this%irpt(ip)
    particle%ip = ip
    particle%name = this%name(ip)
    particle%istopweaksink = this%istopweaksink(ip)
    particle%istopzone = this%istopzone(ip)
    particle%idrymeth = this%idrymeth(ip)
    particle%icp = 0
    particle%icu = this%icu(ip)
    particle%ilay = this%ilay(ip)
    particle%izone = this%izone(ip)
    particle%izp = this%izp(ip)
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
    particle%ifrctrn = this%ifrctrn(ip)
    particle%iexmeth = this%iexmeth(ip)
    particle%extol = this%extol(ip)
    particle%iextend = this%extend(ip)
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
    this%izp(ip) = particle%izp
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
    this%ifrctrn(ip) = particle%ifrctrn
    this%iexmeth(ip) = particle%iexmeth
    this%extol(ip) = particle%extol
    this%extend(ip) = particle%iextend
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

end module ParticleModule

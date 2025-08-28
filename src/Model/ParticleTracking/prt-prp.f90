module PrtPrpModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DEM1, DEM5, DONE, LENFTYPE, LINELENGTH, &
                             LENBOUNDNAME, LENPAKLOC, TABLEFT, TABCENTER, &
                             MNORMAL, DSAME, DEP3, DEP9
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use PrtFmiModule, only: PrtFmiType
  use ParticleModule, only: ParticleType, ParticleStoreType, &
                            create_particle, create_particle_store
  use SimModule, only: count_errors, store_error, store_warning, &
                       store_error_filename
  use SimVariablesModule, only: errmsg, warnmsg
  use ParticleTracksModule, only: ParticleTracksType, &
                                  TRACKHEADER, TRACKDTYPES
  use GeomUtilModule, only: point_in_polygon, get_ijk, get_jk
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_reallocate
  use ParticleReleaseScheduleModule, only: ParticleReleaseScheduleType, &
                                           create_release_schedule
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use ErrorUtilModule, only: pstop
  use MathUtilModule, only: arange, is_close
  use MethodModule, only: LEVEL_MODEL, LEVEL_FEATURE, LEVEL_SUBFEATURE

  implicit none

  private
  public :: PrtPrpType
  public :: prp_create

  character(len=LENFTYPE) :: ftype = 'PRP'
  character(len=16) :: text = '             PRP'

  !> @brief Particle release point (PRP) package
  type, extends(BndExtType) :: PrtPrpType
    ! options
    logical(LGP), pointer :: extend => null() !< extend tracking beyond simulation's end
    logical(LGP), pointer :: frctrn => null() !< force ternary solution for quad grids
    logical(LGP), pointer :: drape => null() !< whether to drape particle to topmost active cell
    logical(LGP), pointer :: localz => null() !< compute z coordinates local to the release cell
    integer(I4B), pointer :: istopweaksink => null() !< weak sink option: 0 = no stop, 1 = stop
    integer(I4B), pointer :: istopzone => null() !< optional stop zone number: 0 = no stop zone
    integer(I4B), pointer :: idrymeth => null() !< dry tracking method: 0 = drop, 1 = stop, 2 = stay
    integer(I4B), pointer :: itrkout => null() !< binary track file
    integer(I4B), pointer :: itrkhdr => null() !< track header file
    integer(I4B), pointer :: itrkcsv => null() !< CSV track file
    integer(I4B), pointer :: irlstls => null() !< release time file
    integer(I4B), pointer :: iexmeth => null() !< method for iterative solution of particle exit location and time in generalized Pollock's method
    integer(I4B), pointer :: ichkmeth => null() !< method for checking particle release coordinates are in the specified cells, 0 = none, 1 = eager
    integer(I4B), pointer :: icycwin => null() !< cycle detection window size
    real(DP), pointer :: extol => null() !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    real(DP), pointer :: rttol => null() !< tolerance for coincident particle release times
    real(DP), pointer :: rtfreq => null() !< frequency for regularly spaced release times
    real(DP), pointer :: offset => null() !< release time offset
    real(DP), pointer :: stoptime => null() !< stop time for all release points
    real(DP), pointer :: stoptraveltime => null() !< stop travel time for all points
    !
    type(PrtFmiType), pointer :: fmi => null() !< flow model interface
    type(ParticleStoreType), pointer :: particles => null() !< particle store
    type(ParticleReleaseScheduleType), pointer :: schedule => null() !< particle release schedule
    integer(I4B), pointer :: nreleasepoints => null() !< number of release points
    integer(I4B), pointer :: nreleasetimes => null() !< number of user-specified particle release times
    integer(I4B), pointer :: nparticles => null() !< number of particles released
    integer(I4B), pointer, contiguous :: rptnode(:) => null() !< release point reduced nns
    integer(I4B), pointer, contiguous :: rptzone(:) => null() !< release point zone numbers
    real(DP), pointer, contiguous :: rptx(:) => null() !< release point x coordinates
    real(DP), pointer, contiguous :: rpty(:) => null() !< release point y coordinates
    real(DP), pointer, contiguous :: rptz(:) => null() !< release point z coordinates
    real(DP), pointer, contiguous :: rptm(:) => null() !< total mass released from point
    character(len=LENBOUNDNAME), pointer, contiguous :: rptname(:) => null() !< release point names
  contains
    procedure :: prp_allocate_arrays
    procedure :: prp_allocate_scalars
    procedure :: bnd_ar => prp_ar
    procedure :: bnd_ad => prp_ad
    procedure :: bnd_rp => prp_rp
    procedure :: bnd_cq_simrate => prp_cq_simrate
    procedure :: bnd_da => prp_da
    procedure :: define_listlabel
    procedure :: prp_set_pointers
    procedure :: source_options => prp_options
    procedure :: source_dimensions => prp_dimensions
    procedure :: prp_log_options
    procedure :: prp_packagedata
    procedure :: prp_releasetimes
    procedure :: prp_load_releasetimefrequency
    procedure :: release
    procedure :: log_release
    procedure :: validate_release_point
    procedure :: initialize_particle
    procedure, public :: bnd_obs_supported => prp_obs_supported
    procedure, public :: bnd_df_obs => prp_df_obs
  end type PrtPrpType

contains

  !> @brief Create a new particle release point package
  subroutine prp_create(packobj, id, ibcnum, inunit, iout, namemodel, &
                        pakname, input_mempath, fmi)
    ! dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: input_mempath
    type(PrtFmiType), pointer :: fmi
    ! local
    type(PrtPrpType), pointer :: prpobj
    ! formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'PRP PARTICLE RELEASE POINT PACKAGE', &
       &' INPUT READ FROM MEMPATH: ', a, /)"

    ! allocate the object and assign values to object variables
    allocate (prpobj)
    packobj => prpobj

    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, input_mempath)
    prpobj%text = text

    ! allocate scalars
    call prpobj%prp_allocate_scalars()

    ! initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 4
    packobj%iscloc = 1

    ! store pointer to flow model interface
    prpobj%fmi => fmi

    ! if prp is enabled, print a message identifying it
    if (inunit > 0) write (iout, fmtheader) input_mempath
  end subroutine prp_create

  !> @brief Deallocate memory
  subroutine prp_da(this)
    class(PrtPrpType) :: this

    ! Deallocate parent
    call this%BndExtType%bnd_da()

    ! Deallocate scalars
    call mem_deallocate(this%localz)
    call mem_deallocate(this%extend)
    call mem_deallocate(this%offset)
    call mem_deallocate(this%stoptime)
    call mem_deallocate(this%stoptraveltime)
    call mem_deallocate(this%istopweaksink)
    call mem_deallocate(this%istopzone)
    call mem_deallocate(this%drape)
    call mem_deallocate(this%idrymeth)
    call mem_deallocate(this%nreleasepoints)
    call mem_deallocate(this%nreleasetimes)
    call mem_deallocate(this%nparticles)
    call mem_deallocate(this%itrkout)
    call mem_deallocate(this%itrkhdr)
    call mem_deallocate(this%itrkcsv)
    call mem_deallocate(this%irlstls)
    call mem_deallocate(this%frctrn)
    call mem_deallocate(this%iexmeth)
    call mem_deallocate(this%ichkmeth)
    call mem_deallocate(this%icycwin)
    call mem_deallocate(this%extol)
    call mem_deallocate(this%rttol)
    call mem_deallocate(this%rtfreq)

    ! Deallocate arrays
    call mem_deallocate(this%rptx)
    call mem_deallocate(this%rpty)
    call mem_deallocate(this%rptz)
    call mem_deallocate(this%rptnode)
    call mem_deallocate(this%rptm)
    call mem_deallocate(this%rptname, 'RPTNAME', this%memoryPath)

    ! Deallocate objects
    call this%particles%destroy(this%memoryPath)
    call this%schedule%destroy()
    deallocate (this%particles)
    deallocate (this%schedule)
  end subroutine prp_da

  !> @ brief Set pointers to model variables
  subroutine prp_set_pointers(this, ibound, izone)
    class(PrtPrpType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    integer(I4B), dimension(:), pointer, contiguous :: izone

    this%ibound => ibound
    this%rptzone => izone
  end subroutine prp_set_pointers

  !> @brief Allocate arrays
  subroutine prp_allocate_arrays(this, nodelist, auxvar)
    ! dummy
    class(PrtPrpType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! local
    integer(I4B) :: nps

    call this%BndExtType%allocate_arrays()

    ! Allocate particle store, starting with the number
    ! of release points (arrays resized if/when needed)
    call create_particle_store( &
      this%particles, &
      this%nreleasepoints, &
      this%memoryPath)

    ! Allocate arrays
    call mem_allocate(this%rptx, this%nreleasepoints, 'RPTX', this%memoryPath)
    call mem_allocate(this%rpty, this%nreleasepoints, 'RPTY', this%memoryPath)
    call mem_allocate(this%rptz, this%nreleasepoints, 'RPTZ', this%memoryPath)
    call mem_allocate(this%rptm, this%nreleasepoints, 'RPTMASS', &
                      this%memoryPath)
    call mem_allocate(this%rptnode, this%nreleasepoints, 'RPTNODER', &
                      this%memoryPath)
    call mem_allocate(this%rptname, LENBOUNDNAME, this%nreleasepoints, &
                      'RPTNAME', this%memoryPath)

    ! Initialize arrays
    do nps = 1, this%nreleasepoints
      this%rptm(nps) = DZERO
    end do
  end subroutine prp_allocate_arrays

  !> @brief Allocate scalars
  subroutine prp_allocate_scalars(this)
    class(PrtPrpType) :: this

    ! Allocate parent's scalars
    call this%BndExtType%allocate_scalars()

    ! Allocate scalars for this type
    call mem_allocate(this%localz, 'LOCALZ', this%memoryPath)
    call mem_allocate(this%extend, 'EXTEND', this%memoryPath)
    call mem_allocate(this%offset, 'OFFSET', this%memoryPath)
    call mem_allocate(this%stoptime, 'STOPTIME', this%memoryPath)
    call mem_allocate(this%stoptraveltime, 'STOPTRAVELTIME', this%memoryPath)
    call mem_allocate(this%istopweaksink, 'ISTOPWEAKSINK', this%memoryPath)
    call mem_allocate(this%istopzone, 'ISTOPZONE', this%memoryPath)
    call mem_allocate(this%drape, 'DRAPE', this%memoryPath)
    call mem_allocate(this%idrymeth, 'IDRYMETH', this%memoryPath)
    call mem_allocate(this%nreleasepoints, 'NRELEASEPOINTS', this%memoryPath)
    call mem_allocate(this%nreleasetimes, 'NRELEASETIMES', this%memoryPath)
    call mem_allocate(this%nparticles, 'NPARTICLES', this%memoryPath)
    call mem_allocate(this%itrkout, 'ITRKOUT', this%memoryPath)
    call mem_allocate(this%itrkhdr, 'ITRKHDR', this%memoryPath)
    call mem_allocate(this%itrkcsv, 'ITRKCSV', this%memoryPath)
    call mem_allocate(this%irlstls, 'IRLSTLS', this%memoryPath)
    call mem_allocate(this%frctrn, 'FRCTRN', this%memoryPath)
    call mem_allocate(this%iexmeth, 'IEXMETH', this%memoryPath)
    call mem_allocate(this%ichkmeth, 'ICHKMETH', this%memoryPath)
    call mem_allocate(this%icycwin, 'ICYCWIN', this%memoryPath)
    call mem_allocate(this%extol, 'EXTOL', this%memoryPath)
    call mem_allocate(this%rttol, 'RTTOL', this%memoryPath)
    call mem_allocate(this%rtfreq, 'RTFREQ', this%memoryPath)

    ! Set values
    this%localz = .false.
    this%extend = .false.
    this%offset = DZERO
    this%stoptime = huge(1d0)
    this%stoptraveltime = huge(1d0)
    this%istopweaksink = 0
    this%istopzone = 0
    this%drape = .false.
    this%idrymeth = 0
    this%nreleasepoints = 0
    this%nreleasetimes = 0
    this%nparticles = 0
    this%itrkout = 0
    this%itrkhdr = 0
    this%itrkcsv = 0
    this%irlstls = 0
    this%frctrn = .false.
    this%iexmeth = 0
    this%ichkmeth = 1
    this%icycwin = 0
    this%extol = DEM5
    this%rttol = DSAME * DEP9
    this%rtfreq = DZERO

  end subroutine prp_allocate_scalars

  !> @ brief Allocate and read period data
  subroutine prp_ar(this)
    ! dummy variables
    class(PrtPrpType), intent(inout) :: this
    ! local variables
    integer(I4B) :: n

    call this%obs%obs_ar()

    if (this%inamedbound /= 0) then
      do n = 1, this%nreleasepoints
        this%boundname(n) = this%rptname(n)
      end do
    end if
    do n = 1, this%nreleasepoints
      this%nodelist(n) = this%rptnode(n)
    end do
  end subroutine prp_ar

  !> @brief Advance a time step and release particles if scheduled.
  subroutine prp_ad(this)
    use TdisModule, only: totalsimtime
    class(PrtPrpType) :: this
    integer(I4B) :: ip, it
    real(DP) :: t

    ! Notes
    ! -----
    ! Each release point can be thought of as
    ! a gumball machine with infinite supply:
    ! a point can release an arbitrary number
    ! of particles, but only one at any time.
    ! Coincident release times are merged to
    ! a single time by the release scheduler.

    ! Reset mass accumulators for this time step.
    do ip = 1, this%nreleasepoints
      this%rptm(ip) = DZERO
    end do

    ! Advance the release schedule and check if
    ! any releases will be made this time step.
    call this%schedule%advance()
    if (.not. this%schedule%any()) return

    ! Log the schedule to the list file.
    call this%log_release()

    ! Expand the particle store. We know from the
    ! schedule how many particles will be released.
    call this%particles%resize( &
      this%particles%num_stored() + &
      (this%nreleasepoints * this%schedule%count()), &
      this%memoryPath)

    ! Release a particle from each point for
    ! each release time in the current step.
    do ip = 1, this%nreleasepoints
      do it = 1, this%schedule%count()
        t = this%schedule%times(it)
        ! Skip the release time if it's before the simulation
        ! starts, or if no `extend_tracking`, after it ends.
        if (t < DZERO) then
          write (warnmsg, '(a,g0,a)') &
            'Skipping negative release time (t=', t, ').'
          call store_warning(warnmsg)
          cycle
        else if (t > totalsimtime .and. .not. this%extend) then
          write (warnmsg, '(a,g0,a)') &
            'Skipping release time falling after the end of the &
            &simulation (t=', t, '). Enable EXTEND_TRACKING to &
            &release particles after the simulation end time.'
          call store_warning(warnmsg)
          cycle
        end if
        call this%release(ip, t)
      end do
    end do
  end subroutine prp_ad

  !> @brief Log the release scheduled for this time step.
  subroutine log_release(this)
    class(PrtPrpType), intent(inout) :: this !< prp
    if (this%iprpak > 0) then
      write (this%iout, "(1x,/1x,a,1x,i0)") &
        'PARTICLE RELEASE FOR PRP', this%ibcnum
      call this%schedule%log(this%iout)
    end if
  end subroutine log_release

  !> @brief Verify that the release point is in the cell.
  !!
  !! Terminate with an error if the release point lies outside the
  !! given cell, or if the point is above or below the grid top or
  !! bottom, respectively.
  !<
  subroutine validate_release_point(this, ic, x, y, z)
    class(PrtPrpType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: ic !< cell index
    real(DP), intent(in) :: x, y, z !< release point
    ! local
    real(DP), allocatable :: polyverts(:, :)

    call this%fmi%dis%get_polyverts(ic, polyverts)
    if (.not. point_in_polygon(x, y, polyverts)) then
      write (errmsg, '(a,g0,a,g0,a,i0)') &
        'Error: release point (x=', x, ', y=', y, ') is not in cell ', &
        this%dis%get_nodeuser(ic)
      call store_error(errmsg, terminate=.false.)
      call store_error_filename(this%input_fname)
    end if
    if (z > maxval(this%dis%top)) then
      write (errmsg, '(a,g0,a,g0,a,i0)') &
        'Error: release point (z=', z, ') is above grid top ', &
        maxval(this%dis%top)
      call store_error(errmsg, terminate=.false.)
      call store_error_filename(this%input_fname)
    else if (z < minval(this%dis%bot)) then
      write (errmsg, '(a,g0,a,g0,a,i0)') &
        'Error: release point (z=', z, ') is below grid bottom ', &
        minval(this%dis%bot)
      call store_error(errmsg, terminate=.false.)
      call store_error_filename(this%input_fname)
    end if
    deallocate (polyverts)
  end subroutine validate_release_point

  !> Release a particle at the specified time.
  !!
  !! Releasing a particle entails validating the particle's
  !! coordinates and settings, transforming its coordinates
  !! if needed, initializing the particle's initial tracking
  !! time to the given release time, storing the particle in
  !! the particle store (from which the PRT model will later
  !! retrieve it, apply the tracking method, and check it in
  !! again), and accumulating the particle's mass (the total
  !! mass released from each release point is calculated for
  !! budget reporting).
  !<
  subroutine release(this, ip, trelease)
    ! dummy
    class(PrtPrpType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: ip !< particle index
    real(DP), intent(in) :: trelease !< release time
    ! local
    integer(I4B) :: np
    type(ParticleType), pointer :: particle

    call this%initialize_particle(particle, ip, trelease)
    np = this%nparticles + 1
    this%nparticles = np
    call this%particles%put(particle, np)
    deallocate (particle)
    this%rptm(ip) = this%rptm(ip) + DONE ! TODO configurable mass

  end subroutine release

  subroutine initialize_particle(this, particle, ip, trelease)
    use ParticleModule, only: TERM_UNRELEASED
    class(PrtPrpType), intent(inout) :: this !< this instance
    type(ParticleType), pointer, intent(inout) :: particle !< the particle
    integer(I4B), intent(in) :: ip !< particle index
    real(DP), intent(in) :: trelease !< release time
    ! local
    integer(I4B) :: irow, icol, ilay, icpl
    integer(I4B) :: ic, icu, ic_old
    real(DP) :: x, y, z
    real(DP) :: top, bot, hds

    ic = this%rptnode(ip)
    icu = this%dis%get_nodeuser(ic)

    call create_particle(particle)

    if (size(this%boundname) /= 0) then
      particle%name = this%boundname(ip)
    else
      particle%name = ''
    end if

    particle%irpt = ip
    particle%istopweaksink = this%istopweaksink
    particle%istopzone = this%istopzone
    particle%idrymeth = this%idrymeth
    particle%icu = icu

    select type (dis => this%dis)
    type is (DisType)
      call get_ijk(icu, dis%nrow, dis%ncol, dis%nlay, irow, icol, ilay)
    type is (DisvType)
      call get_jk(icu, dis%ncpl, dis%nlay, icpl, ilay)
    end select
    particle%ilay = ilay
    particle%izone = this%rptzone(ic)
    particle%istatus = 0 ! status 0 until tracking starts
    ! If the cell is inactive, either drape the particle
    ! to the top-most active cell beneath it if drape is
    ! enabled, or else terminate permanently unreleased.
    if (this%ibound(ic) == 0) then
      ic_old = ic
      if (this%drape) then
        call this%dis%highest_active(ic, this%ibound)
        if (ic == ic_old .or. this%ibound(ic) == 0) then
          ! negative unreleased status signals to the
          ! tracking method that we haven't yet saved
          ! a termination record, it needs to do so.
          particle%istatus = -1 * TERM_UNRELEASED
        end if
      else
        particle%istatus = -1 * TERM_UNRELEASED
      end if
    end if

    ! Load coordinates and transform if needed
    x = this%rptx(ip)
    y = this%rpty(ip)
    if (this%localz) then
      top = this%fmi%dis%top(ic)
      bot = this%fmi%dis%bot(ic)
      hds = this%fmi%gwfhead(ic)
      z = bot + this%rptz(ip) * (hds - bot)
    else
      z = this%rptz(ip)
    end if

    if (this%ichkmeth > 0) &
      call this%validate_release_point(ic, x, y, z)

    particle%x = x
    particle%y = y
    particle%z = z
    particle%trelease = trelease

    ! Set stop time to earlier of STOPTIME and STOPTRAVELTIME
    if (this%stoptraveltime == huge(1d0)) then
      particle%tstop = this%stoptime
    else
      particle%tstop = particle%trelease + this%stoptraveltime
      if (this%stoptime < particle%tstop) particle%tstop = this%stoptime
    end if

    particle%ttrack = particle%trelease
    particle%itrdomain(LEVEL_MODEL) = 0
    particle%iboundary(LEVEL_MODEL) = 0
    particle%itrdomain(LEVEL_FEATURE) = ic
    particle%iboundary(LEVEL_FEATURE) = 0
    particle%itrdomain(LEVEL_SUBFEATURE) = 0
    particle%iboundary(LEVEL_SUBFEATURE) = 0
    particle%frctrn = this%frctrn
    particle%iexmeth = this%iexmeth
    particle%extend = this%extend
    particle%icycwin = this%icycwin
    particle%extol = this%extol
  end subroutine initialize_particle

  !> @ brief Read and prepare period data for particle input
  subroutine prp_rp(this)
    ! modules
    use TdisModule, only: kper, nper
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    ! dummy variables
    class(PrtPrpType), intent(inout) :: this
    ! local variables
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: settings
    integer(I4B), pointer :: iper, ionper, nlist
    character(len=LINELENGTH), allocatable :: lines(:)
    integer(I4B) :: n

    ! set pointer to last and next period loaded
    call mem_setptr(iper, 'IPER', this%input_mempath)
    call mem_setptr(ionper, 'IONPER', this%input_mempath)

    if (kper == 1 .and. &
        (iper == 0) .and. &
        (ionper > nper) .and. &
        size(this%schedule%time_select%times) == 0) then
      ! If the user hasn't provided any release settings (neither
      ! explicit release times, release time frequency, or period
      ! block release settings), default to a single release at the
      ! start of the first period's first time step.
      allocate (lines(1))
      lines(1) = "FIRST"
      call this%schedule%advance(lines=lines)
      deallocate (lines)
      return
    else if (iper /= kper) then
      return
    end if

    ! set input context pointers
    call mem_setptr(nlist, 'NBOUND', this%input_mempath)
    call mem_setptr(settings, 'SETTING', this%input_mempath)

    ! allocate and set input
    allocate (lines(nlist))
    do n = 1, nlist
      lines(n) = settings(n)
    end do

    ! update schedule
    if (size(lines) > 0) &
      call this%schedule%advance(lines=lines)

    ! cleanup
    deallocate (lines)
  end subroutine prp_rp

  !> @ brief Calculate flow between package and model.
  subroutine prp_cq_simrate(this, hnew, flowja, imover)
    ! modules
    use TdisModule, only: delt
    ! dummy variables
    class(PrtPrpType) :: this
    real(DP), dimension(:), intent(in) :: hnew
    real(DP), dimension(:), intent(inout) :: flowja !< flow between package and model
    integer(I4B), intent(in) :: imover !< flag indicating if the mover package is active
    ! local variables
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: idiag
    real(DP) :: rrate

    ! If no boundaries, skip flow calculations.
    if (this%nbound <= 0) return

    ! Loop through each boundary calculating flow.
    do i = 1, this%nbound
      node = this%nodelist(i)
      rrate = DZERO
      ! If cell is no-flow or constant-head, then ignore it.
      if (node > 0) then
        ! Calculate the flow rate into the cell.
        idiag = this%dis%con%ia(node)
        rrate = this%rptm(i) * (DONE / delt) ! reciprocal of tstp length
        flowja(idiag) = flowja(idiag) + rrate
      end if

      ! Save simulated value to simvals array.
      this%simvals(i) = rrate
    end do
  end subroutine prp_cq_simrate

  subroutine define_listlabel(this)
    class(PrtPrpType), intent(inout) :: this
    ! not implemented, not used
  end subroutine define_listlabel

  !> @brief Indicates whether observations are supported.
  logical function prp_obs_supported(this)
    class(PrtPrpType) :: this
    prp_obs_supported = .true.
  end function prp_obs_supported

  !> @brief Store supported observations
  subroutine prp_df_obs(this)
    ! dummy
    class(PrtPrpType) :: this
    ! local
    integer(I4B) :: indx
    call this%obs%StoreObsType('prp', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor

    ! Store obs type and assign procedure pointer
    ! for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine prp_df_obs

  !> @ brief Set options specific to PrtPrpType
  subroutine prp_options(this)
    ! -- modules
    use ConstantsModule, only: LENVARNAME, DZERO, MNORMAL
    use MemoryManagerExtModule, only: mem_set_value
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: getunit, openfile
    use PrtPrpInputModule, only: PrtPrpParamFoundType
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this
    ! -- local variables
    character(len=LENVARNAME), dimension(3) :: drytrack_method = &
      &[character(len=LENVARNAME) :: 'DROP', 'STOP', 'STAY']
    character(len=LENVARNAME), dimension(2) :: coorcheck_method = &
      &[character(len=LENVARNAME) :: 'NONE', 'EAGER']
    character(len=LINELENGTH) :: trackfile, trackcsvfile, fname
    type(PrtPrpParamFoundType) :: found

    ! -- source base class options
    call this%BndExtType%source_options()

    ! -- update defaults from input context
    call mem_set_value(this%stoptime, 'STOPTIME', this%input_mempath, &
                       found%stoptime)
    call mem_set_value(this%stoptraveltime, 'STOPTRAVELTIME', &
                       this%input_mempath, found%stoptraveltime)
    call mem_set_value(this%istopweaksink, 'ISTOPWEAKSINK', this%input_mempath, &
                       found%istopweaksink)
    call mem_set_value(this%istopzone, 'ISTOPZONE', this%input_mempath, &
                       found%istopzone)
    call mem_set_value(this%drape, 'DRAPE', this%input_mempath, &
                       found%drape)
    call mem_set_value(this%idrymeth, 'IDRYMETH', this%input_mempath, &
                       drytrack_method, found%idrymeth)
    call mem_set_value(trackfile, 'TRACKFILE', this%input_mempath, &
                       found%trackfile)
    call mem_set_value(trackcsvfile, 'TRACKCSVFILE', this%input_mempath, &
                       found%trackcsvfile)
    call mem_set_value(this%localz, 'LOCALZ', this%input_mempath, &
                       found%localz)
    call mem_set_value(this%extend, 'EXTEND', this%input_mempath, &
                       found%extend)
    call mem_set_value(this%extol, 'EXTOL', this%input_mempath, &
                       found%extol)
    call mem_set_value(this%rttol, 'RTTOL', this%input_mempath, &
                       found%rttol)
    call mem_set_value(this%rtfreq, 'RTFREQ', this%input_mempath, &
                       found%rtfreq)
    call mem_set_value(this%frctrn, 'FRCTRN', this%input_mempath, &
                       found%frctrn)
    call mem_set_value(this%iexmeth, 'IEXMETH', this%input_mempath, &
                       found%iexmeth)
    call mem_set_value(this%ichkmeth, 'ICHKMETH', this%input_mempath, &
                       coorcheck_method, found%ichkmeth)
    call mem_set_value(this%icycwin, 'ICYCWIN', this%input_mempath, found%icycwin)

    ! update internal state and validate input
    if (found%idrymeth) then
      if (this%idrymeth == 0) then
        write (errmsg, '(a)') 'Unsupported dry tracking method. &
          &DRY_TRACKING_METHOD must be "DROP", "STOP", or "STAY"'
        call store_error(errmsg)
      else
        ! adjust for method zero indexing
        this%idrymeth = this%idrymeth - 1
      end if
    end if

    if (found%extol) then
      if (this%extol <= DZERO) &
        call store_error('EXIT_SOLVE_TOLERANCE MUST BE POSITIVE')
    end if

    if (found%rttol) then
      if (this%rttol <= DZERO) &
        call store_error('RELEASE_TIME_TOLERANCE MUST BE POSITIVE')
    end if

    if (found%rtfreq) then
      if (this%rtfreq <= DZERO) &
        call store_error('RELEASE_TIME_FREQUENCY MUST BE POSITIVE')
    end if

    if (found%iexmeth) then
      if (.not. (this%iexmeth /= 1 .or. this%iexmeth /= 2)) &
        call store_error('DEV_EXIT_SOLVE_METHOD MUST BE &
          &1 (BRENT) OR 2 (CHANDRUPATLA)')
    end if

    if (found%ichkmeth) then
      if (this%ichkmeth == 0) then
        write (errmsg, '(a)') 'Unsupported coordinate check method. &
          &COORDINATE_CHECK_METHOD must be "NONE" or "EAGER"'
        call store_error(errmsg)
      else
        ! adjust for method zero based indexing
        this%ichkmeth = this%ichkmeth - 1
      end if
    end if

    if (found%icycwin) then
      if (this%icycwin < 0) &
        call store_error('CYCLE_DETECTION_WINDOW MUST BE NON-NEGATIVE')
    end if

    ! fileout options
    if (found%trackfile) then
      this%itrkout = getunit()
      call openfile(this%itrkout, this%iout, trackfile, 'DATA(BINARY)', &
                    form, access, filstat_opt='REPLACE', &
                    mode_opt=MNORMAL)
      ! open and write ascii header spec file
      this%itrkhdr = getunit()
      fname = trim(trackfile)//'.hdr'
      call openfile(this%itrkhdr, this%iout, fname, 'CSV', &
                    filstat_opt='REPLACE', mode_opt=MNORMAL)
      write (this%itrkhdr, '(a,/,a)') TRACKHEADER, TRACKDTYPES
    end if

    if (found%trackcsvfile) then
      this%itrkcsv = getunit()
      call openfile(this%itrkcsv, this%iout, trackcsvfile, 'CSV', &
                    filstat_opt='REPLACE')
      write (this%itrkcsv, '(a)') TRACKHEADER
    end if

    ! terminate if any errors were detected
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! log found options
    call this%prp_log_options(found, trackfile, trackcsvfile)

    ! Create release schedule now that we know
    ! the coincident release time tolerance
    this%schedule => create_release_schedule(tolerance=this%rttol)
  end subroutine prp_options

  !> @ brief Log options specific to PrtPrpType
  subroutine prp_log_options(this, found, trackfile, trackcsvfile)
    ! -- modules
    use PrtPrpInputModule, only: PrtPrpParamFoundType
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this
    type(PrtPrpParamFoundType), intent(in) :: found
    character(len=*), intent(in) :: trackfile
    character(len=*), intent(in) :: trackcsvfile
    ! -- local variables
    ! formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    write (this%iout, '(1x,a)') 'PROCESSING PARTICLE INPUT DIMENSIONS'

    if (found%frctrn) then
      write (this%iout, '(4x,a)') &
        'IF DISV, TRACKING WILL USE THE TERNARY METHOD REGARDLESS OF CELL TYPE'
    end if

    if (found%trackfile) then
      write (this%iout, fmttrkbin) trim(adjustl(trackfile)), this%itrkout
    end if

    if (found%trackcsvfile) then
      write (this%iout, fmttrkcsv) trim(adjustl(trackcsvfile)), this%itrkcsv
    end if

    write (this%iout, '(1x,a)') 'END OF PARTICLE INPUT DIMENSIONS'
  end subroutine prp_log_options

  !> @ brief Set dimensions specific to PrtPrpType
  subroutine prp_dimensions(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use PrtPrpInputModule, only: PrtPrpParamFoundType
    ! -- dummy variables
    class(PrtPrpType), intent(inout) :: this
    ! -- local variables
    type(PrtPrpParamFoundType) :: found

    call mem_set_value(this%nreleasepoints, 'NRELEASEPTS', this%input_mempath, &
                       found%nreleasepts)
    call mem_set_value(this%nreleasetimes, 'NRELEASETIMES', this%input_mempath, &
                       found%nreleasetimes)

    write (this%iout, '(1x,a)') 'PROCESSING PARTICLE INPUT DIMENSIONS'
    write (this%iout, '(4x,a,i0)') 'NRELEASEPTS = ', this%nreleasepoints
    write (this%iout, '(4x,a,i0)') 'NRELEASETIMES = ', this%nreleasetimes
    write (this%iout, '(1x,a)') 'END OF PARTICLE INPUT DIMENSIONS'

    ! set maxbound and nbound to nreleasepts
    this%maxbound = this%nreleasepoints
    this%nbound = this%nreleasepoints

    ! allocate arrays for prp package
    call this%prp_allocate_arrays()

    ! read packagedata and releasetimes blocks
    call this%prp_packagedata()
    call this%prp_releasetimes()
    call this%prp_load_releasetimefrequency()
  end subroutine prp_dimensions

  !> @brief Load package data (release points).
  subroutine prp_packagedata(this)
    use MemoryManagerModule, only: mem_setptr
    use GeomUtilModule, only: get_node
    use CharacterStringModule, only: CharacterStringType
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    integer(I4B), dimension(:), pointer, contiguous :: irptno
    integer(I4B), dimension(:, :), pointer, contiguous :: cellids
    real(DP), dimension(:), pointer, contiguous :: xrpts, yrpts, zrpts
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: boundnames
    character(len=LENBOUNDNAME) :: bndName, bndNameTemp
    character(len=9) :: cno
    integer(I4B), dimension(:), allocatable :: nboundchk
    integer(I4B), dimension(:), pointer :: cellid
    integer(I4B) :: n, noder, nodeu, rptno

    ! set input context pointers
    call mem_setptr(irptno, 'IRPTNO', this%input_mempath)
    call mem_setptr(cellids, 'CELLID', this%input_mempath)
    call mem_setptr(xrpts, 'XRPT', this%input_mempath)
    call mem_setptr(yrpts, 'YRPT', this%input_mempath)
    call mem_setptr(zrpts, 'ZRPT', this%input_mempath)
    call mem_setptr(boundnames, 'BOUNDNAME', this%input_mempath)

    ! allocate and initialize temporary variables
    allocate (nboundchk(this%nreleasepoints))
    do n = 1, this%nreleasepoints
      nboundchk(n) = 0
    end do

    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%packName)) &
      //' PACKAGEDATA'

    do n = 1, size(irptno)

      rptno = irptno(n)

      if (rptno < 1 .or. rptno > this%nreleasepoints) then
        write (errmsg, '(a,i0,a,i0,a)') &
          'Expected ', this%nreleasepoints, ' release points. &
          &Points must be numbered from 1 to ', this%nreleasepoints, '.'
        call store_error(errmsg)
        cycle
      end if

      ! increment nboundchk
      nboundchk(rptno) = nboundchk(rptno) + 1

      ! set cellid
      cellid => cellids(:, n)

      ! set node user
      if (this%dis%ndim == 1) then
        nodeu = cellid(1)
      elseif (this%dis%ndim == 2) then
        nodeu = get_node(cellid(1), 1, cellid(2), &
                         this%dis%mshape(1), 1, &
                         this%dis%mshape(2))
      else
        nodeu = get_node(cellid(1), cellid(2), cellid(3), &
                         this%dis%mshape(1), &
                         this%dis%mshape(2), &
                         this%dis%mshape(3))
      end if

      ! set noder
      noder = this%dis%get_nodenumber(nodeu, 1)
      if (noder <= 0) then
        cycle
      else
        this%rptnode(rptno) = noder
      end if

      if (this%localz .and. (zrpts(n) < 0 .or. zrpts(n) > 1)) then
        call store_error('Local z coordinate must fall in the interval [0, 1]')
        cycle
      end if

      ! set coordinates
      this%rptx(rptno) = xrpts(n)
      this%rpty(rptno) = yrpts(n)
      this%rptz(rptno) = zrpts(n)

      ! set default boundname
      write (cno, '(i9.9)') rptno
      bndName = 'PRP'//cno

      ! read boundnames from file, if provided
      if (this%inamedbound /= 0) then
        bndNameTemp = boundnames(n)
        if (bndNameTemp /= '') bndName = bndNameTemp
      else
        bndName = ''
      end if

      ! set boundname
      this%rptname(rptno) = bndName
    end do

    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%packName))//' PACKAGEDATA'

    ! check for duplicate or missing particle release points
    do n = 1, this%nreleasepoints
      if (nboundchk(n) == 0) then
        write (errmsg, '(a,a,1x,i0,a)') 'No data specified for particle ', &
          'release point', n, '.'
        call store_error(errmsg)
      else if (nboundchk(n) > 1) then
        write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
          'Data for particle release point', n, 'specified', nboundchk(n), &
          'times.'
        call store_error(errmsg)
      end if
    end do

    ! terminate if any errors were detected
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! cleanup
    deallocate (nboundchk)
  end subroutine prp_packagedata

  !> @brief Load explicitly specified release times.
  subroutine prp_releasetimes(this)
    use MemoryManagerModule, only: mem_setptr, get_isize
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    real(DP), dimension(:), pointer, contiguous :: time
    integer(I4B) :: n, isize
    real(DP), allocatable :: times(:)

    if (this%nreleasetimes <= 0) return

    ! allocate times array
    allocate (times(this%nreleasetimes))

    ! check if input array was read
    call get_isize('TIME', this%input_mempath, isize)

    if (isize <= 0) then
      errmsg = "RELEASTIMES block expected when &
        &NRELEASETIMES dimension is non-zero."
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    ! set input context pointer
    call mem_setptr(time, 'TIME', this%input_mempath)

    ! set input data
    do n = 1, size(time)
      times(n) = time(n)
    end do

    ! register times with the release schedule
    call this%schedule%time_select%extend(times)

    ! make sure times strictly increase
    if (.not. this%schedule%time_select%increasing()) then
      errmsg = "RELEASTIMES block entries must strictly increase."
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    ! deallocate
    deallocate (times)
  end subroutine prp_releasetimes

  !> @brief Load regularly spaced release times if configured.
  subroutine prp_load_releasetimefrequency(this)
    ! modules
    use TdisModule, only: totalsimtime
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    real(DP), allocatable :: times(:)

    ! check if a release time frequency is configured
    if (this%rtfreq <= DZERO) return

    ! create array of regularly-spaced release times
    times = arange( &
            start=DZERO, &
            stop=totalsimtime, &
            step=this%rtfreq)

    ! register times with release schedule
    call this%schedule%time_select%extend(times)

    ! make sure times strictly increase
    if (.not. this%schedule%time_select%increasing()) then
      errmsg = "Release times must strictly increase"
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    ! deallocate
    deallocate (times)

  end subroutine prp_load_releasetimefrequency

end module PrtPrpModule

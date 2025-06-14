module PrtPrpModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DEM1, DEM5, DONE, LENFTYPE, LINELENGTH, &
                             LENBOUNDNAME, LENPAKLOC, TABLEFT, TABCENTER, &
                             MNORMAL, DSAME, DEP3, DEP9
  use BndModule, only: BndType
  use ObsModule, only: DefaultObsIdProcessor
  use TableModule, only: TableType, table_cr
  use TimeSeriesModule, only: TimeSeriesType
  use TimeSeriesRecordModule, only: TimeSeriesRecordType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  GetTimeSeriesLinkFromList
  use BlockParserModule, only: BlockParserType
  use PrtFmiModule, only: PrtFmiType
  use ParticleModule, only: ParticleType, ParticleStoreType, &
                            create_particle, create_particle_store
  use SimModule, only: count_errors, store_error, store_error_unit, &
                       store_warning
  use SimVariablesModule, only: errmsg, warnmsg
  use TrackControlModule, only: TrackControlType
  use GeomUtilModule, only: point_in_polygon, get_ijk, get_jk
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_reallocate
  use ReleaseScheduleModule, only: ReleaseScheduleType, create_release_schedule
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use ErrorUtilModule, only: pstop
  use MathUtilModule, only: arange, is_close
  use ArrayHandlersModule, only: ExpandArray

  implicit none

  private
  public :: PrtPrpType
  public :: prp_create

  character(len=LENFTYPE) :: ftype = 'PRP'
  character(len=16) :: text = '             PRP'

  !> @brief Particle release point (PRP) package
  type, extends(BndType) :: PrtPrpType
    type(PrtFmiType), pointer :: fmi => null() !< flow model interface
    type(ParticleStoreType), pointer :: particles => null() !< particle store
    type(TrackControlType), pointer :: trackctl => null() !< track control
    type(ReleaseScheduleType), pointer :: schedule !< particle release schedule
    integer(I4B), pointer :: nreleasepoints => null() !< number of release points
    integer(I4B), pointer :: nreleasetimes => null() !< number of user-specified particle release times
    integer(I4B), pointer :: nparticles => null() !< number of particles released
    integer(I4B), pointer :: istopweaksink => null() !< weak sink option: 0 = no stop, 1 = stop
    integer(I4B), pointer :: istopzone => null() !< optional stop zone number: 0 = no stop zone
    integer(I4B), pointer :: idrape => null() !< drape option: 0 = do not drape, 1 = drape to topmost active cell
    integer(I4B), pointer :: idrymeth => null() !< dry tracking method: 0 = drop, 1 = stop, 2 = stay
    integer(I4B), pointer :: itrkout => null() !< binary track file
    integer(I4B), pointer :: itrkhdr => null() !< track header file
    integer(I4B), pointer :: itrkcsv => null() !< CSV track file
    integer(I4B), pointer :: irlstls => null() !< release time file
    integer(I4B), pointer :: ilocalz => null() !< compute z coordinates local to the cell
    integer(I4B), pointer :: iextend => null() !< extend tracking beyond simulation's end
    integer(I4B), pointer :: ifrctrn => null() !< force ternary solution for quad grids
    integer(I4B), pointer :: iexmeth => null() !< method for iterative solution of particle exit location and time in generalized Pollock's method
    real(DP), pointer :: extol => null() !< tolerance for iterative solution of particle exit location and time in generalized Pollock's method
    real(DP), pointer :: rttol => null() !< tolerance for coincident particle release times
    real(DP), pointer :: rtfreq => null() !< frequency for regularly spaced release times
    real(DP), pointer :: offset => null() !< release time offset
    real(DP), pointer :: stoptime => null() !< stop time for all release points
    real(DP), pointer :: stoptraveltime => null() !< stop travel time for all points
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
    procedure :: bnd_options => prp_options
    procedure :: read_dimensions => prp_read_dimensions
    procedure :: prp_read_packagedata
    procedure :: prp_read_releasetimes
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
                        pakname, fmi)
    ! dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    type(PrtFmiType), pointer :: fmi
    ! local
    type(PrtPrpType), pointer :: prpobj
    ! formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'PRP PARTICLE RELEASE POINT PACKAGE', &
       &' INPUT READ FROM UNIT ', i0, /)"

    ! allocate the object and assign values to object variables
    allocate (prpobj)
    packobj => prpobj

    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
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
    if (inunit > 0) write (iout, fmtheader) inunit
  end subroutine prp_create

  !> @brief Deallocate memory
  subroutine prp_da(this)
    class(PrtPrpType) :: this

    ! Deallocate parent
    call this%BndType%bnd_da()

    ! Deallocate scalars
    call mem_deallocate(this%ilocalz)
    call mem_deallocate(this%iextend)
    call mem_deallocate(this%offset)
    call mem_deallocate(this%stoptime)
    call mem_deallocate(this%stoptraveltime)
    call mem_deallocate(this%istopweaksink)
    call mem_deallocate(this%istopzone)
    call mem_deallocate(this%idrape)
    call mem_deallocate(this%idrymeth)
    call mem_deallocate(this%nreleasepoints)
    call mem_deallocate(this%nreleasetimes)
    call mem_deallocate(this%nparticles)
    call mem_deallocate(this%itrkout)
    call mem_deallocate(this%itrkhdr)
    call mem_deallocate(this%itrkcsv)
    call mem_deallocate(this%irlstls)
    call mem_deallocate(this%ifrctrn)
    call mem_deallocate(this%iexmeth)
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
    call this%schedule%deallocate()
    deallocate (this%particles)
    deallocate (this%schedule)
  end subroutine prp_da

  !> @ brief Set pointers to model variables
  subroutine prp_set_pointers(this, ibound, izone, trackctl)
    class(PrtPrpType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    integer(I4B), dimension(:), pointer, contiguous :: izone
    type(TrackControlType), pointer :: trackctl

    this%ibound => ibound
    this%rptzone => izone
    this%trackctl => trackctl
  end subroutine prp_set_pointers

  !> @brief Allocate arrays
  subroutine prp_allocate_arrays(this, nodelist, auxvar)
    ! dummy
    class(PrtPrpType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! local
    integer(I4B) :: nps

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
    call mem_allocate(this%rptm, this%nreleasepoints, 'RPTMASS', this%memoryPath)
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
    call this%BndType%allocate_scalars()

    ! Allocate scalars for this type
    call mem_allocate(this%ilocalz, 'ILOCALZ', this%memoryPath)
    call mem_allocate(this%iextend, 'IEXTEND', this%memoryPath)
    call mem_allocate(this%offset, 'OFFSET', this%memoryPath)
    call mem_allocate(this%stoptime, 'STOPTIME', this%memoryPath)
    call mem_allocate(this%stoptraveltime, 'STOPTRAVELTIME', this%memoryPath)
    call mem_allocate(this%istopweaksink, 'ISTOPWEAKSINK', this%memoryPath)
    call mem_allocate(this%istopzone, 'ISTOPZONE', this%memoryPath)
    call mem_allocate(this%idrape, 'IDRAPE', this%memoryPath)
    call mem_allocate(this%idrymeth, 'IDRYMETH', this%memoryPath)
    call mem_allocate(this%nreleasepoints, 'NRELEASEPOINTS', this%memoryPath)
    call mem_allocate(this%nreleasetimes, 'NRELEASETIMES', this%memoryPath)
    call mem_allocate(this%nparticles, 'NPARTICLES', this%memoryPath)
    call mem_allocate(this%itrkout, 'ITRKOUT', this%memoryPath)
    call mem_allocate(this%itrkhdr, 'ITRKHDR', this%memoryPath)
    call mem_allocate(this%itrkcsv, 'ITRKCSV', this%memoryPath)
    call mem_allocate(this%irlstls, 'IRLSTLS', this%memoryPath)
    call mem_allocate(this%ifrctrn, 'IFRCTRN', this%memoryPath)
    call mem_allocate(this%iexmeth, 'IEXMETH', this%memoryPath)
    call mem_allocate(this%extol, 'EXTOL', this%memoryPath)
    call mem_allocate(this%rttol, 'RTTOL', this%memoryPath)
    call mem_allocate(this%rtfreq, 'RTFREQ', this%memoryPath)

    ! Set values
    this%ilocalz = 0
    this%iextend = 0
    this%offset = DZERO
    this%stoptime = huge(1d0)
    this%stoptraveltime = huge(1d0)
    this%istopweaksink = 0
    this%istopzone = 0
    this%idrape = 0
    this%idrymeth = 0
    this%nreleasepoints = 0
    this%nreleasetimes = 0
    this%nparticles = 0
    this%itrkout = 0
    this%itrkhdr = 0
    this%itrkcsv = 0
    this%irlstls = 0
    this%ifrctrn = 0
    this%iexmeth = 0
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
    call this%BndType%allocate_arrays()
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
        else if (t > totalsimtime .and. this%iextend == 0) then
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
      call store_error_unit(this%inunit, terminate=.true.)
    end if
    if (z > maxval(this%dis%top)) then
      write (errmsg, '(a,g0,a,g0,a,i0)') &
        'Error: release point (z=', z, ') is above grid top ', &
        maxval(this%dis%top)
      call store_error(errmsg, terminate=.false.)
      call store_error_unit(this%inunit, terminate=.true.)
    else if (z < minval(this%dis%bot)) then
      write (errmsg, '(a,g0,a,g0,a,i0)') &
        'Error: release point (z=', z, ') is below grid bottom ', &
        minval(this%dis%bot)
      call store_error(errmsg, terminate=.false.)
      call store_error_unit(this%inunit, terminate=.true.)
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
    use ParticleModule, only: TERM_UNRELEASED, LVL_MODEL, LVL_CELL, LVL_SUBCELL
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
      if (this%idrape > 0) then
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
    if (this%ilocalz > 0) then
      top = this%fmi%dis%top(ic)
      bot = this%fmi%dis%bot(ic)
      hds = this%fmi%gwfhead(ic)
      z = bot + this%rptz(ip) * (hds - bot)
    else
      z = this%rptz(ip)
    end if

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
    particle%itrdomain(LVL_MODEL) = 0
    particle%iboundary(LVL_MODEL) = 0
    particle%itrdomain(LVL_CELL) = ic
    particle%iboundary(LVL_CELL) = 0
    particle%itrdomain(LVL_SUBCELL) = 0
    particle%iboundary(LVL_SUBCELL) = 0
    particle%ifrctrn = this%ifrctrn
    particle%iexmeth = this%iexmeth
    particle%iextend = this%iextend
    particle%extol = this%extol
  end subroutine initialize_particle

  !> @ brief Read and prepare period data for particle input
  subroutine prp_rp(this)
    ! modules
    use TdisModule, only: kper, nper
    use InputOutputModule, only: urword
    ! dummy variables
    class(PrtPrpType), intent(inout) :: this
    ! local variables
    integer(I4B) :: ierr
    logical(LGP) :: is_found
    logical(LGP) :: end_of_block
    logical(LGP) :: no_blocks
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH), allocatable :: lines(:)
    ! formats
    character(len=*), parameter :: fmtblkerr = &
                      "('Looking for BEGIN PERIOD iper.  &
                      &Found ', a, ' instead.')"
    character(len=*), parameter :: fmt_steps = &
                                   "(6x,'TIME STEP(S) ',50(I0,' '))" ! 50 limit is similar to STEPS in OC
    character(len=*), parameter :: fmt_freq = &
                                   "(6x,'EVERY ',I0,' TIME STEP(S)')"
    character(len=*), parameter :: fmt_fracs = &
                                   "(6x,50(f10.3,' '))"

    ! Set ionper to the stress period number for which a new block of data
    ! will be read.
    if (this%inunit == 0) return

    ! get stress period data
    no_blocks = .false.
    if (this%ionper < kper) then
      ! get period block
      call this%parser%GetBlock('PERIOD', is_found, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (is_found) then
        ! read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else
        ! PERIOD block not found
        if (ierr < 0) then
          if (kper == 1) then
            ! End of file found; no period data for the simulation.
            no_blocks = .true.
          else
            ! End of file found; no more period data.
            this%ionper = nper + 1
          end if
        else
          ! Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
    end if

    ! If the user hasn't provided any release settings (neither
    ! explicit release times, release time frequency, or period
    ! block release settings), default to a single release at the
    ! start of the first period's first time step.
    if (no_blocks .and. &
        kper == 1 .and. &
        size(this%schedule%time_select%times) == 0) then
      allocate (lines(1))
      line = "FIRST"
      lines(1) = line
      call this%schedule%advance(lines=lines)
    else if (this%ionper == kper) then
      ! If the current stress period matches the
      ! block we are reading, parse the setting
      ! and register it with the schedule.
      allocate (lines(0))
      recordloop: do
        call this%parser%GetNextLine(end_of_block)
        if (end_of_block) exit recordloop
        call this%parser%GetCurrentLine(line)
        call ExpandArray(lines)
        lines(size(lines)) = line
      end do recordloop
      if (size(lines) > 0) &
        call this%schedule%advance(lines=lines)
      deallocate (lines)
    end if

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

  !> @brief Set options specific to PrtPrpType
  subroutine prp_options(this, option, found)
    use OpenSpecModule, only: access, form
    use ConstantsModule, only: MAXCHARLEN, DZERO
    use InputOutputModule, only: urword, getunit, openfile
    use TrackFileModule, only: TRACKHEADER, TRACKDTYPES
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical(LGP), intent(inout) :: found
    ! locals
    character(len=MAXCHARLEN) :: fname
    character(len=MAXCHARLEN) :: keyword
    ! formats
    character(len=*), parameter :: fmttrkbin = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO BINARY FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"
    character(len=*), parameter :: fmttrkcsv = &
      "(4x, 'PARTICLE TRACKS WILL BE SAVED TO CSV FILE: ', a, /4x, &
    &'OPENED ON UNIT: ', I0)"

    select case (option)
    case ('STOPTIME')
      this%stoptime = this%parser%GetDouble()
      found = .true.
    case ('STOPTRAVELTIME')
      this%stoptraveltime = this%parser%GetDouble()
      found = .true.
    case ('STOP_AT_WEAK_SINK')
      this%istopweaksink = 1
      found = .true.
    case ('ISTOPZONE')
      this%istopzone = this%parser%GetInteger()
      found = .true.
    case ('DRAPE')
      this%idrape = 1
      found = .true.
    case ('DRY_TRACKING_METHOD')
      call this%parser%GetStringCaps(keyword)
      select case (keyword)
      case ('DROP')
        this%idrymeth = 0
      case ('STOP')
        this%idrymeth = 1
      case ('STAY')
        this%idrymeth = 2
      case default
        write (errmsg, '(a, a)') &
          'Unknown dry tracking method: ', trim(keyword)
        call store_error(errmsg)
        write (errmsg, '(a, a)') &
          'DRY must be "DROP", "STOP" or "STAY"'
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end select
      found = .true.
    case ('TRACK')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        ! parse filename
        call this%parser%GetString(fname)
        ! open binary output file
        this%itrkout = getunit()
        call openfile(this%itrkout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, filstat_opt='REPLACE', &
                      mode_opt=MNORMAL)
        write (this%iout, fmttrkbin) trim(adjustl(fname)), this%itrkout
        ! open and write ascii header spec file
        this%itrkhdr = getunit()
        fname = trim(fname)//'.hdr'
        call openfile(this%itrkhdr, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE', mode_opt=MNORMAL)
        write (this%itrkhdr, '(a,/,a)') TRACKHEADER, TRACKDTYPES
      else
        call store_error('OPTIONAL TRACK KEYWORD MUST BE '// &
                         'FOLLOWED BY FILEOUT')
      end if
      found = .true.
    case ('TRACKCSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        ! parse filename
        call this%parser%GetString(fname)
        ! open CSV output file and write headers
        this%itrkcsv = getunit()
        call openfile(this%itrkcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE')
        write (this%iout, fmttrkcsv) trim(adjustl(fname)), this%itrkcsv
        write (this%itrkcsv, '(a)') TRACKHEADER
      else
        call store_error('OPTIONAL TRACKCSV KEYWORD MUST BE &
          &FOLLOWED BY FILEOUT')
      end if
      found = .true.
    case ('LOCAL_Z')
      this%ilocalz = 1
      found = .true.
    case ('EXTEND_TRACKING')
      this%iextend = 1
      found = .true.
    case ('EXIT_SOLVE_TOLERANCE')
      this%extol = this%parser%GetDouble()
      if (this%extol <= DZERO) &
        call store_error('EXIT_SOLVE_TOLERANCE MUST BE POSITIVE')
      found = .true.
    case ('RELEASE_TIME_TOLERANCE')
      this%rttol = this%parser%GetDouble()
      if (this%rttol <= DZERO) &
        call store_error('RELEASE_TIME_TOLERANCE MUST BE POSITIVE')
      found = .true.
    case ('RELEASE_TIME_FREQUENCY')
      this%rtfreq = this%parser%GetDouble()
      if (this%rtfreq <= DZERO) &
        call store_error('RELEASE_TIME_FREQUENCY MUST BE POSITIVE')
      found = .true.
    case ('DEV_FORCETERNARY')
      call this%parser%DevOpt()
      this%ifrctrn = 1
      write (this%iout, '(4x,a)') &
        'TRACKING WILL BE DONE USING THE TERNARY METHOD REGARDLESS OF CELL TYPE'
      found = .true.
    case ('DEV_EXIT_SOLVE_METHOD')
      call this%parser%DevOpt()
      this%iexmeth = this%parser%GetInteger()
      if (.not. (this%iexmeth /= 1 .or. this%iexmeth /= 2)) &
        call store_error('DEV_EXIT_SOLVE_METHOD MUST BE &
          &1 (BRENT) OR 2 (CHANDRUPATLA)')
      found = .true.
    case default
      found = .false.
    end select

    ! Catch unrecognized options
    if (.not. found) then
      errmsg = "UNKNOWN PRP OPTION '"//trim(keyword)//"'."
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if

    ! Create release schedule now that we know
    ! the coincident release time tolerance
    this%schedule => create_release_schedule(tol=this%rttol)

  end subroutine prp_options

  !> @brief Read package dimensions
  subroutine prp_read_dimensions(this)
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock

    ! get dimension block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)

    ! parse dimension block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING PARTICLE INPUT DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('NRELEASEPTS')
          this%nreleasepoints = this%parser%GetInteger()
        case ('NRELEASETIMES')
          this%nreleasetimes = this%parser%GetInteger()
        case default
          write (errmsg, &
                 '(4x,a,a)') '****ERROR. UNKNOWN PARTICLE INPUT DIMENSION: ', &
            trim(keyword)
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF PARTICLE INPUT DIMENSIONS'
    else
      call store_error('ERROR.  REQUIRED DIMENSIONS BLOCK NOT FOUND.')
    end if

    ! set maxbound and nbound to nreleasepts
    this%maxbound = this%nreleasepoints
    this%nbound = this%nreleasepoints

    ! allocate arrays for prp package
    call this%prp_allocate_arrays()

    ! read packagedata and releasetimes blocks
    call this%prp_read_packagedata()
    call this%prp_read_releasetimes()
    call this%prp_load_releasetimefrequency()
  end subroutine prp_read_dimensions

  !> @brief Load package data (release points).
  subroutine prp_read_packagedata(this)
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    character(len=LINELENGTH) :: cellid
    character(len=LENBOUNDNAME) :: bndName, bndNameTemp
    character(len=9) :: cno
    logical :: isfound
    logical :: endOfBlock
    integer(I4B) :: ival
    integer(I4B) :: n
    integer(I4B) :: ierr
    character(len=LENBOUNDNAME), dimension(:), allocatable :: nametxt
    integer(I4B), dimension(:), allocatable :: nboundchk
    integer(I4B), dimension(:), allocatable :: noder
    real(DP), dimension(:), allocatable :: x
    real(DP), dimension(:), allocatable :: y
    real(DP), dimension(:), allocatable :: z
    real(DP), dimension(:), allocatable :: tstop
    ! format
    character(len=*), parameter :: fmttend = &
      "('end time (', G0, ') must be greater than or equal to the              &
     &begin time (', G0, ').')"

    ! allocate temporary variables
    allocate (noder(this%nreleasepoints))
    allocate (x(this%nreleasepoints))
    allocate (y(this%nreleasepoints))
    allocate (z(this%nreleasepoints))
    allocate (tstop(this%nreleasepoints))
    allocate (nametxt(this%nreleasepoints))
    allocate (nboundchk(this%nreleasepoints))

    ! initialize temporary variables
    do n = 1, this%nreleasepoints
      nboundchk(n) = 0
    end do

    ! read particle release point data
    ! get particle release points block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportopenclose=.true.)

    ! parse block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%packName)) &
        //' PACKAGEDATA'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        ival = this%parser%GetInteger()
        n = ival

        if (n < 1 .or. n > this%nreleasepoints) then
          write (errmsg, '(a,i0,a,i0,a)') &
            'Expected ', this%nreleasepoints, ' release points. &
            &Points must be numbered from 1 to ', this%nreleasepoints, '.'
          call store_error(errmsg)
          cycle
        end if

        ! increment nboundchk
        nboundchk(n) = nboundchk(n) + 1

        ! node number
        call this%parser%GetCellid(this%dis%ndim, cellid)
        noder(n) = this%dis%noder_from_cellid(cellid, this%inunit, this%iout)

        ! x, y, z coordinates
        x(n) = this%parser%GetDouble()
        y(n) = this%parser%GetDouble()
        z(n) = this%parser%GetDouble()

        if (this%ilocalz > 0 .and. (z(n) < 0 .or. z(n) > 1)) then
          call store_error('Local z coordinate must fall in the interval [0, 1]')
          cycle
        end if

        ! set default boundname
        write (cno, '(i9.9)') n
        bndName = 'PRP'//cno

        ! read boundnames from file, if provided
        if (this%inamedbound /= 0) then
          call this%parser%GetStringCaps(bndNameTemp)
          if (bndNameTemp /= '') &
            bndName = bndNameTemp
        else
          bndName = ''
        end if

        ! store temp boundnames
        nametxt(n) = bndName
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
    else
      call store_error('Required packagedata block not found.')
    end if

    ! terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if

    ! fill particle release point data with data stored in temporary local arrays
    do n = 1, this%nreleasepoints
      this%rptnode(n) = noder(n)
      this%rptx(n) = x(n)
      this%rpty(n) = y(n)
      this%rptz(n) = z(n)
      this%rptname(n) = nametxt(n)
    end do

    ! deallocate local storage
    deallocate (noder)
    deallocate (x)
    deallocate (y)
    deallocate (z)
    deallocate (tstop)
    deallocate (nametxt)
    deallocate (nboundchk)
  end subroutine prp_read_packagedata

  !> @brief Load explicitly specified release times.
  subroutine prp_read_releasetimes(this)
    ! dummy
    class(PrtPrpType), intent(inout) :: this
    ! local
    integer(I4B) :: i, ierr
    logical(LGP) :: eob, found, success
    real(DP) :: t
    real(DP), allocatable :: times(:)

    ! get releasetimes block
    call this%parser%GetBlock('RELEASETIMES', found, ierr, &
                              supportOpenClose=.true., &
                              blockRequired=.false.)

    ! raise an error if releasetimes has a dimension
    ! but no block was found, otherwise return early
    if (.not. found) then
      if (this%nreleasetimes <= 0) return
      write (errmsg, '(a, i0)') &
        "Expected RELEASETIMES with length ", this%nreleasetimes
      call store_error(errmsg)
      call this%parser%StoreErrorUnit(terminate=.true.)
    end if

    ! allocate times array
    allocate (times(this%nreleasetimes))

    ! read times from the block
    write (this%iout, '(/1x,a)') &
      'PROCESSING '//trim(adjustl(this%text))//' RELEASETIMES'
    do i = 1, this%nreleasetimes
      call this%parser%GetNextLine(eob)
      if (eob) exit
      call this%parser%TryGetDouble(t, success)
      if (.not. success) then
        errmsg = "Failed to read double precision value"
        call store_error(errmsg)
        call this%parser%StoreErrorUnit(terminate=.true.)
      end if
      times(i) = t
    end do

    ! register times with the release schedule
    call this%schedule%time_select%extend(times)

    ! make sure times strictly increase
    if (.not. this%schedule%time_select%increasing()) then
      errmsg = "Release times must strictly increase"
      call store_error(errmsg)
      call this%parser%StoreErrorUnit(terminate=.true.)
    end if

    ! deallocate
    deallocate (times)

  end subroutine prp_read_releasetimes

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
      call this%parser%StoreErrorUnit(terminate=.true.)
    end if

    ! deallocate
    deallocate (times)

  end subroutine prp_load_releasetimefrequency

end module PrtPrpModule

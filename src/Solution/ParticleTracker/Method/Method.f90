!> @brief Particle tracking strategies
module MethodModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO
  use ErrorUtilModule, only: pstop
  use SubcellModule, only: SubcellType
  use ParticleModule, only: ParticleType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use ParticleEventModule, only: ParticleEventType
  use ReleaseEventModule, only: ReleaseEventType
  use TimeStepEventModule, only: TimeStepEventType
  use TerminationEventModule, only: TerminationEventType
  use WeakSinkEventModule, only: WeakSinkEventType
  use UserTimeEventModule, only: UserTimeEventType
  use FeatExitEventModule, only: FeatExitEventType
  use BaseDisModule, only: DisBaseType
  use PrtFmiModule, only: PrtFmiType
  use CellModule, only: CellType
  use CellDefnModule, only: CellDefnType
  use TimeSelectModule, only: TimeSelectType
  use MathUtilModule, only: is_close
  use DomainModule, only: DomainType
  use ExitSolutionModule, only: ExitSolutionType
  use ListModule, only: ListType
  implicit none

  public :: LEVEL_MODEL, LEVEL_FEATURE, LEVEL_SUBFEATURE

  !> @brief Tracking method level enumeration.
  !!
  !> Tracking levels: 1: model, 2: grid feature, 3: grid subfeature.
  !! A tracking level identifies the domain through which a tracking
  !! method is responsible for moving a particle. Methods operate on
  !! a particular level and delegate to submethods for levels higher
  !! than (i.e. below the scope of) their own.
  !<
  enum, bind(C)
    enumerator :: LEVEL_MODEL = 1
    enumerator :: LEVEL_FEATURE = 2
    enumerator :: LEVEL_SUBFEATURE = 3
  end enum

  private
  public :: MethodType

  !> @brief Base type for particle tracking methods.
  !!
  !! The PRT tracking algorithm invokes a "tracking method" for each
  !! domain. A domain can be a model, cell in a model, or subcell in
  !! a cell. Tracking proceeds recursively, delegating to a possibly
  !! arbitrary number of subdomains (currently, only the three above
  !! are recognized). A tracking method is responsible for advancing
  !! a particle through a domain, delegating to subdomains as needed
  !! depending on cell geometry (implementing the strategy pattern).
  !<
  type, abstract :: MethodType
    character(len=40), pointer, public :: name !< method name
    logical(LGP), public :: delegates !< whether the method delegates
    type(PrtFmiType), pointer, public :: fmi => null() !< ptr to fmi
    class(CellType), pointer, public :: cell => null() !< ptr to the current cell
    class(SubcellType), pointer, public :: subcell => null() !< ptr to the current subcell
    type(ParticleEventDispatcherType), pointer, public :: events => null() !< ptr to event dispatcher
    type(TimeSelectType), pointer, public :: tracktimes => null() !< ptr to user-defined tracking times
    integer(I4B), dimension(:), pointer, contiguous, public :: izone => null() !< pointer to zone numbers
    real(DP), dimension(:), pointer, contiguous, public :: flowja => null() !< pointer to intercell flows
    real(DP), dimension(:), pointer, contiguous, public :: porosity => null() !< pointer to aquifer porosity
    real(DP), dimension(:), pointer, contiguous, public :: retfactor => null() !< pointer to retardation factor
  contains
    ! Implemented in all subtypes
    procedure(apply), deferred :: apply !< apply the method to the particle
    procedure(assess), deferred :: assess !< assess conditions before tracking
    procedure(deallocate), deferred :: deallocate !< deallocate the method object
    procedure :: get_level !< get the tracking method level
    ! Overridden in subtypes that delegate
    procedure :: pass !< pass the particle to the next subdomain
    procedure :: load !< load the subdomain tracking method
    procedure :: find_exits !< find domain exit solutions
    procedure :: pick_exit
    ! Implemented here
    procedure :: init
    procedure :: track
    procedure :: try_pass
    ! Event firing methods
    procedure :: release
    procedure :: terminate
    procedure :: timestep
    procedure :: weaksink
    procedure :: usertime
  end type MethodType

  abstract interface
    subroutine apply(this, particle, tmax)
      import DP
      import MethodType
      import ParticleType
      class(MethodType), intent(inout) :: this
      type(ParticleType), pointer, intent(inout) :: particle
      real(DP), intent(in) :: tmax
    end subroutine apply
    subroutine assess(this, particle, cell_defn, tmax)
      import DP
      import MethodType
      import ParticleType
      import CellDefnType
      class(MethodType), intent(inout) :: this
      type(ParticleType), pointer, intent(inout) :: particle
      type(CellDefnType), pointer, intent(inout) :: cell_defn
      real(DP), intent(in) :: tmax
    end subroutine assess
    subroutine deallocate (this)
      import MethodType
      class(MethodType), intent(inout) :: this
    end subroutine deallocate
  end interface

contains

  !> @brief Initialize the method with pointers to model data.
  subroutine init(this, fmi, cell, subcell, events, tracktimes, &
                  izone, flowja, porosity, retfactor)
    class(MethodType), intent(inout) :: this
    type(PrtFmiType), intent(in), pointer, optional :: fmi
    class(CellType), intent(in), pointer, optional :: cell
    class(SubcellType), intent(in), pointer, optional :: subcell
    type(ParticleEventDispatcherType), intent(in), pointer, optional :: events
    type(TimeSelectType), intent(in), pointer, optional :: tracktimes
    integer(I4B), intent(in), pointer, optional :: izone(:)
    real(DP), intent(in), pointer, optional :: flowja(:)
    real(DP), intent(in), pointer, optional :: porosity(:)
    real(DP), intent(in), pointer, optional :: retfactor(:)

    if (present(fmi)) this%fmi => fmi
    if (present(cell)) this%cell => cell
    if (present(subcell)) this%subcell => subcell
    if (present(events)) this%events => events
    if (present(tracktimes)) this%tracktimes => tracktimes
    if (present(izone)) this%izone => izone
    if (present(flowja)) this%flowja => flowja
    if (present(porosity)) this%porosity => porosity
    if (present(retfactor)) this%retfactor => retfactor
  end subroutine init

  !> @brief Track the particle over subdomains of the given
  ! level until the particle terminates or tmax is reached.
  recursive subroutine track(this, particle, level, tmax)
    ! dummy
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B) :: level
    real(DP), intent(in) :: tmax
    ! local
    logical(LGP) :: advancing
    integer(I4B) :: nextlevel
    class(methodType), pointer :: submethod

    advancing = .true.
    nextlevel = level + 1
    do while (advancing)
      call this%load(particle, nextlevel, submethod)
      call submethod%apply(particle, tmax)
      call this%try_pass(particle, nextlevel, advancing)
    end do
  end subroutine track

  !> @brief Try passing the particle to the next subdomain.
  subroutine try_pass(this, particle, nextlevel, advancing)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B) :: nextlevel
    logical(LGP) :: advancing

    if (particle%advancing) then
      ! if still advancing, pass to the next subdomain.
      ! if that puts us on a boundary, then we're done.
      call this%pass(particle)
      if (particle%iboundary(nextlevel - 1) .ne. 0) &
        advancing = .false.
    else
      ! otherwise we're already done so
      ! reset the domain boundary value.
      advancing = .false.
      particle%iboundary = 0
    end if
  end subroutine try_pass

  !> @brief Get tracking method level.
  function get_level(this) result(level)
    class(MethodType), intent(in) :: this
    integer(I4B) :: level
    level = -1 ! suppress compiler warning
    call pstop(1, "get_level must be overridden")
  end function get_level

  !> @brief Load subdomain tracking method (submethod).
  subroutine load(this, particle, next_level, submethod)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod
    call pstop(1, "load must be overridden")
  end subroutine load

  !> @brief Pass particle to the next subdomain or to a domain boundary.
  subroutine pass(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    call pstop(1, "pass must be overridden")
  end subroutine pass

  !> @brief Compute candidate exit solutions.
  subroutine find_exits(this, particle, domain)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(DomainType), intent(in) :: domain
    if (.not. this%delegates) &
      call pstop(1, "find_exits called on non-delegating method")
    call pstop(1, "find_exits must be overridden in delegating methods")
  end subroutine find_exits

  !> @brief Choose an exit solution among candidates.
  function pick_exit(this, particle) result(exit_soln)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B) :: exit_soln
    exit_soln = 0 ! suppress compiler warning
    if (.not. this%delegates) &
      call pstop(1, "pick_exit called on non-delegating method")
    call pstop(1, "pick_exit must be overridden in delegating methods")
  end function pick_exit

  !> @brief A particle is released.
  subroutine release(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event
    allocate (ReleaseEventType :: event)
    call this%events%dispatch(particle, event)
    deallocate (event)
  end subroutine release

  !> @brief A particle terminates.
  subroutine terminate(this, particle, status)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B), intent(in), optional :: status
    class(ParticleEventType), pointer :: event
    particle%advancing = .false.
    if (present(status)) particle%istatus = status
    allocate (TerminationEventType :: event)
    call this%events%dispatch(particle, event)
    deallocate (event)
  end subroutine terminate

  !> @brief A time step ends.
  subroutine timestep(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event
    allocate (TimeStepEventType :: event)
    call this%events%dispatch(particle, event)
    deallocate (event)
  end subroutine timestep

  !> @brief A particle leaves a weak sink.
  subroutine weaksink(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event
    allocate (WeakSinkEventType :: event)
    call this%events%dispatch(particle, event)
    deallocate (event)
  end subroutine weaksink

  !> @brief A user-defined tracking time occurs.
  subroutine usertime(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event
    allocate (UserTimeEventType :: event)
    call this%events%dispatch(particle, event)
    deallocate (event)
  end subroutine usertime

end module MethodModule

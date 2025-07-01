!> @brief Particle tracking strategies
module MethodModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO
  use ErrorUtilModule, only: pstop
  use SubcellModule, only: SubcellType
  use ParticleModule, only: ParticleType
  use ParticleEventsModule, only: ParticleEventDispatcherType
  use ParticleEventModule, only: ParticleEventType
  use BaseDisModule, only: DisBaseType
  use PrtFmiModule, only: PrtFmiType
  use CellModule, only: CellType
  use CellDefnModule, only: CellDefnType
  use TimeSelectModule, only: TimeSelectType
  use MathUtilModule, only: is_close
  implicit none

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
    procedure(apply), deferred :: apply
    procedure(deallocate), deferred :: deallocate
    ! Overridden in subtypes that delegate
    procedure :: pass
    procedure :: load
    ! Implemented here
    procedure :: init
    procedure :: dispatch
    procedure :: track
    procedure :: try_pass
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
    subroutine deallocate (this)
      import MethodType
      class(MethodType), intent(inout) :: this
    end subroutine deallocate
  end interface

contains

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

  !> @brief Track the particle over domains of the given
  ! level until the particle terminates or tmax elapses.
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

    ! Advance the particle over subdomains
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

    ! if the particle is done advancing, reset the domain boundary flag.
    if (.not. particle%advancing) then
      particle%iboundary = 0
      advancing = .false.
    else
      ! otherwise pass the particle to the next subdomain.
      ! if that leaves it on a boundary, stop advancing.
      call this%pass(particle)
      if (particle%iboundary(nextlevel - 1) .ne. 0) &
        advancing = .false.
    end if
  end subroutine try_pass

  !> @brief Load the subdomain tracking method (submethod).
  subroutine load(this, particle, next_level, submethod)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod
    call pstop(1, "load must be overridden")
  end subroutine load

  !> @brief Pass the particle to the next subdomain.
  subroutine pass(this, particle)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    call pstop(1, "pass must be overridden")
  end subroutine pass

  !> @brief Dispatch a particle event to handlers.
  subroutine dispatch(this, particle, event_code)
    class(MethodType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer(I4B), intent(in) :: event_code

    select case (event_code)
    case (0)
      call this%events%release(particle)
    case (1)
      call this%events%cellexit(particle)
    case (2)
      call this%events%timestep(particle)
    case (3)
      call this%events%terminate(particle)
    case (4)
      call this%events%weaksink(particle)
    case (5)
      call this%events%usertime(particle)
    case default
      call pstop(1, "Programmer error, bad event code")
    end select
  end subroutine dispatch

end module MethodModule

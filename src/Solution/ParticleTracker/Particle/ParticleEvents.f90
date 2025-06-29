module ParticleEventsModule
  use KindModule, only: DP, I4B, LGP
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
  public :: ExitEventType, TerminationEventType, ReleaseEventType
  public :: TimestepEventType, WeakSinkEventType, UserTimeEventType
  public :: RELEASE, EXIT, TIMESTEP, TERMINATE, WEAKSINK, USERTIME

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

  type, abstract :: ParticleEventType
    integer(I4B) :: code = -1 ! event code above
    integer(I4B) :: kper = 0, kstp = 0
    real(DP) :: time = 0.0_DP
    character(len=40) :: context = ""
    integer(I4B) :: level = 0 ! 1=model, 2=cell, 3=subcell
  contains
    procedure :: get_code
  end type ParticleEventType

  type :: EventWrapperType
    class(ParticleEventType), allocatable :: event
  end type EventWrapperType

  type, extends(ParticleEventType) :: ExitEventType
    integer(I4B), allocatable :: domain(:)
    integer(I4B), allocatable :: boundary(:)
  end type ExitEventType

  type, extends(ParticleEventType) :: TerminationEventType
    character(len=80) :: message = ""
  end type TerminationEventType

  type, extends(ParticleEventType) :: ReleaseEventType
    real(DP) :: coords(3) = 0.0_DP
    integer(I4B) :: domain(3) = 0
  end type ReleaseEventType

  type, extends(ParticleEventType) :: TimestepEventType
  end type TimestepEventType

  type, extends(ParticleEventType) :: WeakSinkEventType
  end type WeakSinkEventType

  type, extends(ParticleEventType) :: UserTimeEventType
  end type UserTimeEventType

contains

  function get_code(this) result(code)
    class(ParticleEventType), intent(in) :: this
    integer(I4B) :: code

    select type (this)
    type is (ReleaseEventType); code = 0
    type is (ExitEventType); code = 1
    type is (TimestepEventType); code = 2
    type is (TerminationEventType); code = 3
    type is (WeakSinkEventType); code = 4
    type is (UserTimeEventType); code = 5
    class default; code = -1
    end select
  end function get_code

end module ParticleEventsModule

module ParticleEventModule
  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: ParticleEventType
  public :: CellExitEventType, TerminationEventType, ReleaseEventType
  public :: TimeStepEventType, WeakSinkEventType, UserTimeEventType
  public :: RELEASE, CELLEXIT, TIMESTEP, TERMINATE, WEAKSINK, USERTIME

  !> @brief Particle event enumeration.
  !!
  !! A number of events may occur to particles, each of which may (or may
  !! not) be of interest to the user. The user selects events to report.
  !<
  enum, bind(C)
    enumerator :: RELEASE = 0 !< particle was released
    enumerator :: CELLEXIT = 1 !< particle exited a cell
    enumerator :: TIMESTEP = 2 !< time step ended
    enumerator :: TERMINATE = 3 !< particle terminated
    enumerator :: WEAKSINK = 4 !< particle exited a weak sink
    enumerator :: USERTIME = 5 !< user-specified tracking time
  end enum

  !> @brief Base type for particle events.
  !!
  !! Events may be identical except for their type/code, reflecting the
  !! fact that several events of interest may occur at a given moment.
  type, abstract :: ParticleEventType
    integer(I4B) :: imdl, iprp, irpt ! release model, package, and point
    real(DP) :: trelease = 0.0_DP ! release time
    integer(I4B) :: kper = 0, kstp = 0 ! period and step
    integer(I4B) :: ilay, icu, izone = 0
    real(DP) :: ttrack = 0.0_DP ! simulation time
    real(DP) :: x = 0.0_DP, y = 0.0_DP, z = 0.0_DP ! particle position
    integer(I4B) :: istatus = -1 ! status code
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: log
  end type ParticleEventType

  type, extends(ParticleEventType) :: CellExitEventType
  end type CellExitEventType

  type, extends(ParticleEventType) :: TerminationEventType
  end type TerminationEventType

  type, extends(ParticleEventType) :: ReleaseEventType
  end type ReleaseEventType

  type, extends(ParticleEventType) :: TimeStepEventType
  end type TimestepEventType

  type, extends(ParticleEventType) :: WeakSinkEventType
  end type WeakSinkEventType

  type, extends(ParticleEventType) :: UserTimeEventType
  end type UserTimeEventType

contains
  integer function get_code(this) result(code)
    class(ParticleEventType), intent(in) :: this

    select type (this)
    type is (ReleaseEventType); code = 0
    type is (CellExitEventType); code = 1
    type is (TimeStepEventType); code = 2
    type is (TerminationEventType); code = 3
    type is (WeakSinkEventType); code = 4
    type is (UserTimeEventType); code = 5
    class default; call pstop(1, "unknown event type")
    end select
  end function get_code

  function get_verb(this) result(str)
    class(ParticleEventType), intent(in) :: this
    character(len=:), allocatable :: str

    select type (this)
    type is (ReleaseEventType); str = "released"
    type is (CellExitEventType); str = "exited cell"
    type is (TimeStepEventType); str = "completed timestep"
    type is (TerminationEventType); str = "terminated"
    type is (WeakSinkEventType); str = "exited weak sink"
    type is (UserTimeEventType); str = "user-specified tracking time"
    class default; call pstop(1, "unknown event type")
    end select
  end function get_verb

  subroutine log(this, iun)
    class(ParticleEventType), intent(inout) :: this
    integer(I4B), intent(in) :: iun

    if (iun >= 0) &
      write (iun, '(*(G0))') &
      'Particle (Model: ', this%imdl, &
      ', Package: ', this%iprp, &
      ', Point: ', this%irpt, &
      ', Time: ', this%trelease, &
      ') ', this%get_verb(), &
      ' in (Layer: ', this%ilay, &
      ', Cell: ', this%icu, &
      ', Zone: ', this%izone, &
      ') at (X: ', this%x, &
      ', Y: ', this%y, &
      ', Z: ', this%z, &
      ', Time: ', this%ttrack, &
      ', Period: ', this%kper, &
      ', Timestep: ', this%kstp, &
      ') with (Status: ', this%istatus, ')'
  end subroutine log

end module ParticleEventModule

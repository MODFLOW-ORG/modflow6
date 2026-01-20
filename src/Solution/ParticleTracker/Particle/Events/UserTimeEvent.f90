module UserTimeEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, USERTIME
  implicit none

  private
  public :: UserTimeEventType

  type, extends(ParticleEventType) :: UserTimeEventType
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_text
  end type UserTimeEventType

contains

  function get_code(this) result(code)
    class(UserTimeEventType), intent(in) :: this
    integer(I4B) :: code
    code = USERTIME
  end function get_code

  function get_verb(this) result(verb)
    class(UserTimeEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'user-specified tracking time'
  end function get_verb

  function get_text(this) result(text)
    class(UserTimeEventType), intent(in) :: this
    character(len=:), allocatable :: text
    character(len=LENHUGELINE) :: temp

    write (temp, '(*(G0))') &
      'Particle from model ', this%imdl, &
      ', package ', this%iprp, &
      ', point ', this%irpt, &
      ', time ', this%trelease, &
      ' '//this%get_verb()// &
      ' in layer ', this%ilay, &
      ', cell ', this%icu, &
      ', zone ', this%izone, &
      ' at x ', this%x, &
      ', y ', this%y, &
      ', z ', this%z, &
      ', time ', this%ttrack, &
      ', period ', this%kper, &
      ', timestep ', this%kstp, &
      ' with status ', this%istatus
    text = trim(adjustl(temp))
  end function get_text

end module UserTimeEventModule

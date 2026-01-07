module WeakSinkEventModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENHUGELINE
  use ErrorUtilModule, only: pstop
  use ParticleModule, only: ParticleType
  use ParticleEventModule, only: ParticleEventType, WEAKSINK
  implicit none

  private
  public :: WeakSinkEventType

  type, extends(ParticleEventType) :: WeakSinkEventType
  contains
    procedure :: get_code
    procedure :: get_verb
    procedure :: get_text
  end type WeakSinkEventType

contains

  function get_code(this) result(code)
    class(WeakSinkEventType), intent(in) :: this
    integer(I4B) :: code
    code = WEAKSINK
  end function get_code

  function get_verb(this) result(verb)
    class(WeakSinkEventType), intent(in) :: this
    character(len=:), allocatable :: verb
    verb = 'exited weak sink'
  end function get_verb

  function get_text(this) result(text)
    class(WeakSinkEventType), intent(in) :: this
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

end module WeakSinkEventModule

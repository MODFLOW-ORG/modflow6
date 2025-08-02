module MethodSubcellModule
  use KindModule, only: DP, I4B
  use MethodModule, only: MethodType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType

  private
  public :: MethodSubcellType

  type, abstract, extends(MethodType) :: MethodSubcellType
  contains
    procedure, public :: assess
  end type MethodSubcellType

contains

  subroutine assess(this, particle, cell_defn, tmax)
    ! dummy
    class(MethodSubcellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! noop
  end subroutine assess

end module MethodSubcellModule
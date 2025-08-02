module MethodModelModule
  use KindModule, only: DP, I4B
  use MethodModule, only: MethodType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType

  private
  public :: MethodModelType

  type, abstract, extends(MethodType) :: MethodModelType
  contains
    procedure, public :: assess
  end type MethodModelType

contains

  subroutine assess(this, particle, cell_defn, tmax)
    ! dummy
    class(MethodModelType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! noop
  end subroutine assess

end module MethodModelModule
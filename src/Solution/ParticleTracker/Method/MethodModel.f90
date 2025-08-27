module MethodModelModule
  use KindModule, only: DP, I4B
  use MethodModule, only: MethodType, LEVEL_MODEL
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType

  private
  public :: MethodModelType

  type, abstract, extends(MethodType) :: MethodModelType
  contains
    procedure, public :: assess
    procedure, public :: get_level
  end type MethodModelType

contains

  !> @brief Check particle reporting/termination status
  subroutine assess(this, particle, cell_defn, tmax)
    ! dummy
    class(MethodModelType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(CellDefnType), pointer, intent(inout) :: cell_defn
    real(DP), intent(in) :: tmax
    ! noop
  end subroutine assess

  !> @brief Get the model method's level.
  function get_level(this) result(level)
    class(MethodModelType), intent(in) :: this
    integer(I4B) :: level
    level = LEVEL_MODEL
  end function get_level

end module MethodModelModule
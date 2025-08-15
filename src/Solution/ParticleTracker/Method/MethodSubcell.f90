module MethodSubcellModule
  use KindModule, only: DP, I4B
  use MethodModule, only: LEVEL_SUBFEATURE
  use MethodCellModule, only: MethodCellType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType
  use SubCellExitEventModule, only: SubCellExitEventType

  private
  public :: MethodSubcellType

  type, abstract, extends(MethodCellType) :: MethodSubcellType
  contains
    procedure, public :: assess
    procedure, public :: subcellexit
    procedure, public :: get_level
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

  !> @brief Particle exits a subcell.
  subroutine subcellexit(this, particle)
    class(MethodSubcellType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    class(ParticleEventType), pointer :: event

    allocate (SubCellExitEventType :: event)
    select type (event)
    type is (SubCellExitEventType)
      event%isc = particle%itrdomain(LEVEL_SUBFEATURE)
      event%exit_face = particle%iboundary(LEVEL_SUBFEATURE)
    end select
    call this%events%dispatch(particle, event)
  end subroutine subcellexit

  !> @brief Get the subcell method's level.
  function get_level(this) result(level)
    class(MethodSubcellType), intent(in) :: this
    integer(I4B) :: level
    level = LEVEL_SUBFEATURE
  end function get_level

end module MethodSubcellModule
module MethodSubcellModule
  use KindModule, only: DP, I4B
  use MethodModule, only: MethodType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType, SubcellExitEventType

  private
  public :: MethodSubcellType

  type, abstract, extends(MethodType) :: MethodSubcellType
  contains
    procedure, public :: assess
    procedure, public :: subcellexit
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

    allocate (SubcellExitEventType :: event)
    select type (event)
    type is (SubcellExitEventType)
      event%isc = particle%idomain(3)
      event%exit_face = particle%iboundary(3)
    end select
    call this%events%dispatch(particle, event)
  end subroutine subcellexit

end module MethodSubcellModule
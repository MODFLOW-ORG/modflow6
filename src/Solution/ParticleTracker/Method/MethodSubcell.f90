module MethodSubcellModule
  use KindModule, only: DP, I4B
  use ErrorUtilModule, only: pstop
  use MethodModule, only: LEVEL_SUBFEATURE
  use MethodCellModule, only: MethodCellType
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType
  use SubCellExitEventModule, only: SubCellExitEventType

  private
  public :: MethodSubcellType

  !> @brief Abstract base type for subcell tracking methods
  type, abstract, extends(MethodCellType) :: MethodSubcellType
  contains
    procedure, public :: assess
    procedure, public :: subcellexit
    procedure, public :: get_level
  end type MethodSubcellType

contains

  !> @brief Assess conditions before tracking
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
    ! local
    integer(I4B) :: i, nhist
    class(*), pointer :: prev

    allocate (SubCellExitEventType :: event)
    select type (event)
    type is (SubCellExitEventType)
      event%isc = particle%itrdomain(LEVEL_SUBFEATURE)
      event%exit_face = particle%iboundary(LEVEL_SUBFEATURE)
    end select
    call this%events%broadcast(particle, event)
    if (particle%icycwin == 0) then
      deallocate (event)
      return
    end if
    if (this%forms_cycle(particle, event)) then
      print *, "Cyclic subcell pathline detected"
      nhist = particle%history%Count()
      do i = 1, nhist
        prev => particle%history%GetItem(i)
        select type (prev)
        class is (ParticleEventType)
          print *, "Back ", nhist - i + 1, ": ", prev%get_text()
        end select
      end do
      print *, "Current :", event%get_text()
      call pstop(1, 'Cyclic subcell pathline detected, aborting')
    else
      call this%store_event(particle, event)
    end if
  end subroutine subcellexit

  !> @brief Get the subcell method's level.
  function get_level(this) result(level)
    class(MethodSubcellType), intent(in) :: this
    integer(I4B) :: level
    level = LEVEL_SUBFEATURE
  end function get_level

end module MethodSubcellModule

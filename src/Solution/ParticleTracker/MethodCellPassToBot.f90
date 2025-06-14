module MethodCellPassToBotModule

  use KindModule, only: DP, I4B
  use MethodCellModule, only: MethodCellType
  use CellDefnModule, only: CellDefnType, create_defn
  use PrtFmiModule, only: PrtFmiType
  use BaseDisModule, only: DisBaseType
  use DisModule, only: DisType
  use DisvModule, only: DisvType
  use ParticleModule, only: ParticleType
  use CellModule, only: CellType
  use SubcellModule, only: SubcellType
  use TrackControlModule, only: TrackControlType
  implicit none

  private
  public :: MethodCellPassToBotType
  public :: create_method_cell_ptb

  type, extends(MethodCellType) :: MethodCellPassToBotType
    private
  contains
    procedure, public :: apply => apply_ptb
    procedure, public :: deallocate
  end type MethodCellPassToBotType

contains

  !> @brief Create a new pass-to-bottom tracking method
  subroutine create_method_cell_ptb(method)
    type(MethodCellPassToBotType), pointer :: method
    allocate (method)
    allocate (method%name)
    method%name = "passtobottom"
    method%delegates = .false.
  end subroutine create_method_cell_ptb

  !> @brief Deallocate the pass-to-bottom tracking method
  subroutine deallocate (this)
    class(MethodCellPassToBotType), intent(inout) :: this
    deallocate (this%name)
  end subroutine deallocate

  !> @brief Pass particle vertically and instantaneously to the cell bottom
  subroutine apply_ptb(this, particle, tmax)
    use ParticleModule, only: TERM_NO_EXITS, LVL_CELL
    ! dummy
    class(MethodCellPassToBotType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    integer(I4B) :: nlay

    ! Check termination/reporting conditions
    call this%check(particle, this%cell%defn, tmax)
    if (.not. particle%advancing) return

    ! Pass to bottom face
    particle%z = this%cell%defn%bot
    particle%iboundary(LVL_CELL) = this%cell%defn%npolyverts + 2

    ! Terminate if in bottom layer
    select type (dis => this%fmi%dis)
    type is (DisType)
      nlay = dis%nlay
    type is (DisvType)
      nlay = dis%nlay
    end select
    if (particle%ilay == nlay) then
      particle%advancing = .false.
      particle%istatus = TERM_NO_EXITS
      call this%save(particle, reason=3)
    end if

    ! Save record
    call this%save(particle, reason=1)
  end subroutine apply_ptb

end module MethodCellPassToBotModule

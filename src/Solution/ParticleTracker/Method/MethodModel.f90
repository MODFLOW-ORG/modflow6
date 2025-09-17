module MethodModelModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE
  use MethodModule, only: MethodType, LEVEL_MODEL
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType
  use ParticleEventModule, only: ParticleEventType
  use MathUtilModule, only: is_Close

  private
  public :: MethodModelType

  type, abstract, extends(MethodType) :: MethodModelType
  contains
    procedure, public :: assess
    procedure, public :: get_level
    ! Utilities
    procedure :: cap_wt_flow
    procedure :: set_no_exit_face
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

  !> @brief Prevent non-boundary upwards flow at the water table.
  !!
  !! Unless the top face is an assigned boundary with outflow,
  !! cells which contain a water table should not have upward
  !! flow through the top (i.e., the water table). Prevent it
  !! by capping the top face flow at zero in these conditions.
  !!
  !! Assumes cell properties and flows are already loaded.
  !<
  subroutine cap_wt_flow(this, defn)
    class(MethodModelType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: defn
    ! local
    integer(I4B) :: itopface

    ! If the cell contains a water table that is not an
    ! assigned boundary face with upward flow, cap the flow
    ! at zero. The cell contains a water table if it's partially
    ! saturated (we know it's not dry because we wouldn't be here)
    ! or if it's saturated and its top neighbor is dry. Saturation
    ! and dryness are determined using threshold-based criteria to
    ! account for possible numerical noise.
    itopface = this%fmi%max_faces ! fmi's lateral face indices are not closed
    if (this%fmi%is_boundary_face(defn%icell, itopface)) return
    if (this%cell_has_water_table(defn)) then
      itopface = defn%npolyverts + 3 ! cell defn's lateral face indices are closed
      defn%faceflow(itopface) = max(DZERO, defn%faceflow(itopface))
    end if

  end subroutine cap_wt_flow

  !> @brief Set flag indicating if the cell has any faces with outflow.
  !! Assumes cell properties and flows are already loaded.
  subroutine set_no_exit_face(this, defn)
    ! dummy
    class(MethodModelType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: defn
    ! local
    integer(I4B) :: m, nfaces

    defn%inoexitface = 1
    nfaces = defn%npolyverts + 3
    do m = 1, nfaces
      if (defn%faceflow(m) < DZERO) defn%inoexitface = 0
    end do

  end subroutine set_no_exit_face

end module MethodModelModule
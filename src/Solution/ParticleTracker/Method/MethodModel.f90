module MethodModelModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE
  use MethodModule, only: MethodType, LEVEL_MODEL
  use ParticleModule, only: ParticleType
  use CellDefnModule, only: CellDefnType, SATURATION_DRY, &
                            SATURATION_WATERTABLE, &
                            SATURATION_SATURATED
  use ParticleEventModule, only: ParticleEventType
  use MathUtilModule, only: is_close
  use ErrorUtilModule, only: pstop

  private
  public :: MethodModelType

  type, abstract, extends(MethodType) :: MethodModelType
  contains
    procedure, public :: assess
    procedure, public :: get_level
    ! cell load utilities
    procedure :: cap_cell_wt_flow
    procedure :: load_cell_no_exit_face
    procedure :: load_cell_saturation_status
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
  subroutine cap_cell_wt_flow(this, defn)
    class(MethodModelType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: defn
    ! local
    integer(I4B) :: itopface

    ! If the cell contains a water table that is not an
    ! assigned boundary face with upflow, cap flow at 0.
    itopface = this%fmi%max_faces ! fmi's lateral face indices are not closed
    if (this%fmi%is_boundary_face(defn%icell, itopface)) return
    if (defn%isatstat == SATURATION_WATERTABLE) then
      itopface = defn%npolyverts + 3 ! cell defn's lateral face indices are closed
      defn%faceflow(itopface) = max(DZERO, defn%faceflow(itopface))
    end if

  end subroutine cap_cell_wt_flow

  !> @brief Set flag indicating if the cell has any faces with outflow.
  !! Assumes cell properties and flows are already loaded.
  subroutine load_cell_no_exit_face(this, defn)
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

  end subroutine load_cell_no_exit_face

  !> @brief Set the saturation status of a cell.
  !! See the status enumeration in CellDefnModule.
  subroutine load_cell_saturation_status(this, defn)
    ! dummy
    class(MethodModelType), intent(inout) :: this
    type(CellDefnType), pointer, intent(inout) :: defn
    ! local
    integer(I4B) :: ic, ictopnbr, itopface, idiag, ipos

    ic = defn%icell
    defn%isatstat = SATURATION_SATURATED

    ! dry?
    if (this%fmi%ibdgwfsat0(ic) == 0) then
      defn%isatstat = SATURATION_DRY
      return
    end if

    ! partially saturated?
    if (this%fmi%gwfsat(ic) < DONE) then
      defn%isatstat = SATURATION_WATERTABLE
      return
    end if

    ! no top neighbor?
    itopface = defn%npolyverts + 3 ! cell defn's lateral face indices are closed
    if (defn%facenbr(itopface) == 0) then
      defn%isatstat = SATURATION_SATURATED
      return
    end if

    ! dry top neighbor?
    idiag = this%fmi%dis%con%ia(ic)
    ipos = idiag + defn%facenbr(itopface)
    ictopnbr = this%fmi%dis%con%ja(ipos)
    if (this%fmi%ibdgwfsat0(ictopnbr) == 0) then
      defn%isatstat = SATURATION_WATERTABLE
      return
    end if

  end subroutine load_cell_saturation_status

end module MethodModelModule

module MethodCellPollockModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO
  use MethodModule, only: MethodType
  use MethodCellModule, only: MethodCellType
  use MethodSubcellPoolModule, only: method_subcell_plck, &
                                     method_subcell_tern
  use CellRectModule, only: CellRectType, create_cell_rect
  use SubcellRectModule, only: SubcellRectType, create_subcell_rect
  use ParticleModule, only: ParticleType
  implicit none

  private
  public :: MethodCellPollockType
  public :: create_method_cell_pollock

  type, extends(MethodCellType) :: MethodCellPollockType
  contains
    procedure, public :: apply => apply_mcp
    procedure, public :: deallocate => destroy_mcp
    procedure, public :: load => load_mcp
    procedure, public :: load_subcell
    procedure, public :: pass => pass_mcp
  end type MethodCellPollockType

contains

  !> @brief Create a tracking method
  subroutine create_method_cell_pollock(method)
    ! dummy
    type(MethodCellPollockType), pointer :: method
    ! local
    type(CellRectType), pointer :: cell
    type(SubcellRectType), pointer :: subcell

    allocate (method)
    call create_cell_rect(cell)
    method%cell => cell
    method%name => method%cell%type
    method%delegates = .true.
    call create_subcell_rect(subcell)
    method%subcell => subcell
  end subroutine create_method_cell_pollock

  !> @brief Destroy the tracking method
  subroutine destroy_mcp(this)
    class(MethodCellPollockType), intent(inout) :: this
    deallocate (this%name)
  end subroutine destroy_mcp

  !> @brief Load subcell tracking method
  subroutine load_mcp(this, particle, next_level, submethod)
    ! modules
    use SubcellModule, only: SubcellType
    ! dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    integer, intent(in) :: next_level
    class(MethodType), pointer, intent(inout) :: submethod

    select type (subcell => this%subcell)
    type is (SubcellRectType)
      call this%load_subcell(particle, subcell)
    end select
    call method_subcell_plck%init( &
      cell=this%cell, &
      subcell=this%subcell, &
      trackctl=this%trackctl, &
      tracktimes=this%tracktimes)
    submethod => method_subcell_plck
    particle%idomain(next_level) = 1
  end subroutine load_mcp

  !> @brief Having exited the lone subcell, pass the particle to the cell face
  !! In this case the lone subcell is the cell.
  subroutine pass_mcp(this, particle)
    ! dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    ! local
    integer(I4B) :: exitface
    integer(I4B) :: entryface

    exitface = particle%iboundary(3)
    ! Map subcell exit face to cell face
    select case (exitface) ! note: exitFace uses Dave's iface convention
    case (0)
      entryface = -1
    case (1)
      entryface = 1
    case (2)
      entryface = 3
    case (3)
      entryface = 4
    case (4)
      entryface = 2
    case (5)
      entryface = 6 ! note: inface=5 same as inface=1 due to wraparound
    case (6)
      entryface = 7
    end select
    if (entryface .eq. -1) then
      particle%iboundary(2) = 0
    else
      if ((entryface .ge. 1) .and. (entryface .le. 4)) then
        ! Account for local cell rotation
        select type (cell => this%cell)
        type is (CellRectType)
          entryface = entryface + cell%ipvOrigin - 1
        end select
        if (entryface .gt. 4) entryface = entryface - 4
      end if
      particle%iboundary(2) = entryface
    end if
  end subroutine pass_mcp

  !> @brief Apply Pollock's method to a rectangular cell
  subroutine apply_mcp(this, particle, tmax)
    ! dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    real(DP), intent(in) :: tmax
    ! local
    real(DP) :: xOrigin
    real(DP) :: yOrigin
    real(DP) :: zOrigin
    real(DP) :: sinrot
    real(DP) :: cosrot

    select type (cell => this%cell)
    type is (CellRectType)
      ! Check termination/reporting conditions
      call this%check(particle, cell%defn, tmax)
      if (.not. particle%advancing) return

      ! Transform model coordinates to local cell coordinates
      ! (translated/rotated but not scaled relative to model)
      xOrigin = cell%xOrigin
      yOrigin = cell%yOrigin
      zOrigin = cell%zOrigin
      sinrot = cell%sinrot
      cosrot = cell%cosrot
      call particle%transform(xOrigin, yOrigin, zOrigin, &
                              sinrot, cosrot)

      ! Track the particle over the cell
      call this%track(particle, 2, tmax)

      ! Transform cell coordinates back to model coordinates
      call particle%transform(xOrigin, yOrigin, zOrigin, &
                              sinrot, cosrot, invert=.true.)
      call particle%reset_transform()
    end select
  end subroutine apply_mcp

  !> @brief Loads the lone rectangular subcell from the rectangular cell
  subroutine load_subcell(this, particle, subcell) !
    ! dummy
    class(MethodCellPollockType), intent(inout) :: this
    type(ParticleType), pointer, intent(inout) :: particle
    type(SubcellRectType), intent(inout) :: subcell

    select type (cell => this%cell)
    type is (CellRectType)
      ! Set subcell number to 1
      subcell%isubcell = 1

      ! Subcell calculations will be done in local subcell coordinates
      subcell%dx = cell%dx
      subcell%dy = cell%dy
      subcell%dz = cell%dz
      subcell%sinrot = DZERO
      subcell%cosrot = DONE
      subcell%xOrigin = DZERO
      subcell%yOrigin = DZERO
      subcell%zOrigin = DZERO

      ! Set subcell edge velocities
      subcell%vx1 = cell%vx1 ! cell velocities already account for retfactor and porosity
      subcell%vx2 = cell%vx2
      subcell%vy1 = cell%vy1
      subcell%vy2 = cell%vy2
      subcell%vz1 = cell%vz1
      subcell%vz2 = cell%vz2
    end select
  end subroutine load_subcell

end module MethodCellPollockModule
